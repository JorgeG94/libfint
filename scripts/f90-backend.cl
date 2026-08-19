;;;; Fortran backend for gen-code.cl.
;;;;
;;;; Load after gen-code.cl, then call gen-f90-cint the way you would gen-cint:
;;;;
;;;;     (load "gen-code.cl")
;;;;     (load "f90-backend.cl")
;;;;     (gen-f90-cint "cint_intor1.f90"
;;;;       '("int1e_kin" (.5 \| p dot p))
;;;;       '("int1e_r"   ( \| r )))
;;;;
;;;; The symbolic layer is untouched.  parser.cl and derivator.cl do not know
;;;; what language they are feeding, and gen-code.cl's combinatorial walk over
;;;; bra and ket operators is reused verbatim -- only the emission differs.
;;;; One change was needed upstream, in gen-code.cl itself: combo-bra,
;;;; combo-ket and combo-opj now go through `emit`, which accepts a function
;;;; where it used to accept only a format string.  See below for why.
;;;;
;;;; THE THREE THINGS THAT DIFFER FROM THE C BACKEND
;;;;
;;;; 1. Pointer arithmetic becomes integer offsets into one flat array.
;;;;    `double *g1 = g0 + envs->g_size*3` is an integer add in both languages;
;;;;    `g1[ix+0]` becomes `g(g1+ix+0)`.  With g declared g(0:) the index
;;;;    arithmetic is character-for-character the C's, which is the point --
;;;;    see PORT_TO_FORTRAN.md 3.6.
;;;;
;;;; 2. The G1E_* macros are not all the same kind of statement.  cpp lets the
;;;;    C write G1E_D_J and G1E_R_J identically even though the first expands
;;;;    to a call and the second to `f = g + envs->g_stride_j`, an assignment.
;;;;    Fortran cannot hide that, so the operator letter has to reach the
;;;;    emitter as a value rather than baked into a format string -- which is
;;;;    what the `emit` change upstream is for.  Getting this wrong produces
;;;;    code that compiles and is wrong, which is how it was found.
;;;;
;;;; 3. No `+=`, so the accumulate form restates its target.
;;;;
;;;; Two-electron integrals are emitted too, including their spinor entry
;;;; points -- gen-f90-int4c2e below.  The four-centre case differs from the
;;;; one-centre one in three ways and no more: there are four centres to walk
;;;; instead of two, the s-accumulation runs over the Rys roots rather than a
;;;; single point, and the driver returns complex for the spinor forms.

;;; ---------------------------------------------------------------- cells

;;;; ------------------------------------------------------------ g pointers
;;;;
;;;; The C's generated gout kernels carry one pointer per g block -- `double
;;;; *g1 = g0 + envs->g_size * 3` -- so a term reads `g1[ix+0]`, one indexed
;;;; load.  Writing that as `g(g1+ix+0)` costs an extra add on every one of
;;;; them, and a gout term is three of them.  Measured on cint_g0_2e_2d,
;;;; which has the same shape by hand, binding the rows took it from 1.39x
;;;; the C to 1.10x.
;;;;
;;;; Bound only when there are few enough blocks to be worth it.  A fourth
;;;; derivative has sixteen; the r-polynomial families have a hundred and
;;;; twenty-eight, and binding those would cost more in descriptors than the
;;;; adds it saves.
(defparameter *f90-gptr-max* 16)
(defvar *f90-gptr* nil)

(defun gref (blk expr)
  (if *f90-gptr*
      (format nil "g~ap(~a)" blk expr)
      (format nil "g(g~a+~a)" blk expr)))

(defun emit-gptr-decls (fout nblk)
  (let ((names (loop for i below nblk collect (format nil "g~ap(:)" i))))
    (loop for chunk on names by (lambda (l) (nthcdr 6 l)) do
          (format fout "      real(dp), pointer, contiguous :: ~{~a~^, ~}~%"
                  (subseq chunk 0 (min 6 (length chunk)))))))

(defun emit-gptr-binds (fout nblk)
  (loop for i below nblk do
        (format fout "      g~ap(0:) => g(g~a:)~%" i i)))

(defun f90-cell-converter (cell fout &optional (with-grids nil))
  "cell-converter with () for [].  Same branching, same order of terms."
  (let ((fac (realpart (phase-of cell)))
        (const@3 (ternary-subscript (consts-of cell)))
        (op@3    (if with-grids
                    (format nil "ig+GRID_BLKSIZE*~a" (ternary-subscript (ops-of cell)))
                    (ternary-subscript (ops-of cell)))))
    (cond ((equal fac 1)
           (cond ((null const@3)
                  (if (null op@3) (format fout " + s(0)") (format fout " + s(~a)" op@3)))
                 ((null op@3) (format fout " + c(~a)*s(0)" const@3))
                 (t (format fout " + c(~a)*s(~a)" const@3 op@3))))
          ((equal fac -1)
           (cond ((null const@3)
                  (if (null op@3) (format fout " - s(0)") (format fout " - s(~a)" op@3)))
                 ((null op@3) (format fout " - c(~a)*s(0)" const@3))
                 (t (format fout " - c(~a)*s(~a)" const@3 op@3))))
          ((< fac 0)
           (cond ((null const@3)
                  (if (null op@3) (format fout " ~a*s(0)" fac) (format fout " ~a*s(~a)" fac op@3)))
                 ((null op@3) (format fout " ~a*c(~a)*s(0)" fac const@3))
                 (t (format fout " ~a*c(~a)*s(~a)" fac const@3 op@3))))
          (t
            (cond ((null const@3)
                   (if (null op@3) (format fout " + ~a*s(0)" fac) (format fout " + ~a*s(~a)" fac op@3)))
                  ((null op@3) (format fout " + ~a*c(~a)*s(0)" fac const@3))
                  (t (format fout " + ~a*c(~a)*s(~a)" fac const@3 op@3)))))))

(defun gen-f90-block (fout flat-script)
  (let ((assemb (to-c-code-string fout #'f90-cell-converter flat-script nil))
        (comp (length flat-script)))
    (loop for s in assemb
          for gid from 0 do
          (emit-gout-line fout (format nil "            gout(n*~a+~a) =" comp gid) s))))

;;; The accumulate form needs the right-hand side PARENTHESISED.
;;;
;;; C's `gout[k] += - s[0] - s[4] - s[8]` evaluates the whole right-hand side
;;; first and then adds it once.  Fortran's `gout(k) = gout(k) - s(0) - s(4)
;;; - s(8)` associates left to right, so it rounds three times against a
;;; running total instead of once against a completed sum.  For a single
;;; primitive the two agree, because the empty branch computes the same
;;; expression; with two or more they differ by an ulp, and int1e_kin -- the
;;; only 1e integral whose gout combines several s terms -- drifted to 1e-13
;;; over a contracted shell pair while every other integral stayed exact.
(defun gen-f90-block+ (fout flat-script)
  (let ((assemb (to-c-code-string fout #'f90-cell-converter flat-script nil))
        (comp (length flat-script)))
    (loop for s in assemb
          for gid from 0 do
          ;; the right-hand side stays parenthesised: C evaluates it as a
          ;; group and adds once, Fortran's restated target would associate
          ;; left to right and round against a running total instead
          (emit-gout-line fout
                          (format nil "            gout(n*~a+~a) = gout(n*~a+~a) +" comp gid comp gid)
                          (format nil " (~a )" s)))))

(defun gen-f90-block-with-empty (fout flat-script)
  (format fout "         if (gout_empty /= 0) then~%")
  (gen-f90-block fout flat-script)
  (format fout "         else~%")
  (gen-f90-block+ fout flat-script)
  (format fout "         end if~%"))

;;; Fortran's 132-column limit is a hard error, not a warning, and four
;;; unrolled Rys roots put this well past it.  Each term goes on its own
;;; continued line: same expression, same left-to-right association, and
;;; legible next to the C it came from.
(defun dump-s-for-nroots-f90 (fout tot-bits nroots &optional (indent "         "))
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "~as(~a) = " indent i)
      (loop for k upto (1- nroots) do
            (when (> k 0) (format fout " &~%~a         " indent))
            (format fout "+ ~a*~a*~a"
                    (gref xbin (format nil "ix+~a" k))
                    (gref ybin (format nil "iy+~a" k))
                    (gref zbin (format nil "iz+~a" k))))
      (format fout "~%"))))

;;; -------------------------------------------------------------- G1E ops

;;; The split the C hides behind cpp.  Calls take (array, offset) pairs; the
;;; shifts are plain offset assignments and cannot be calls at all.
(defparameter *f90-1e-call*
  '(("D_" . "cint_nabla1~a_1e") ("R0" . "cint_x1~a_1e") ("RC" . "cint_x1~a_1e")))
(defparameter *f90-1e-shift*
  '(("R_" . "g_stride_~a")))

;;; Which vector an origin-shifting operation is relative to.
(defparameter *f90-1e-origin*
  '(("R0" . "envs%r~a") ("RC" . "dr~a")))

(defun f90-g1e-emit (fout op tgt src il-expr jl-expr centre)
  (let ((call  (cdr (assoc op *f90-1e-call*  :test #'string=)))
        (shift (cdr (assoc op *f90-1e-shift* :test #'string=)))
        (orig  (cdr (assoc op *f90-1e-origin* :test #'string=)))
        (c (string-downcase centre)))
    (cond
      (call
        (if orig
            ;; One array, two offsets: passing the same array as two dummies,
            ;; one of them written, is illegal aliasing in Fortran even though
            ;; the C spells it as two pointers.
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, 0, envs)~%"
                    (format nil call c) tgt src (format nil orig c) il-expr jl-expr)
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, 0, envs)~%"
                    (format nil call c) tgt src il-expr jl-expr)))
      (shift
        (format fout "      g~a = g~a + envs%~a~%" tgt src (format nil shift c)))
      (t (error "unmapped G1E op ~a" op)))))

;;; Closures matching the argument order the combo walkers use.
(defun f90-fmt-i-fn ()
  (lambda (fout op tgt src lshift)
    (f90-g1e-emit fout op tgt src
                  (format nil "envs%i_l+~a" lshift) "envs%j_l" "i")))

(defun f90-fmt-j-fn ()
  (lambda (fout op tgt src ilen right)
    (f90-g1e-emit fout op tgt src
                  (format nil "envs%i_l+~d" ilen)
                  (format nil "envs%j_l+~a" right) "j")))

(defun f90-fmt-op-fn ()
  (lambda (fout op1 t1 s1 il1 r1 op2 t2 s2 il2 r2 acc-a acc-b)
    (f90-g1e-emit fout op1 t1 s1 (format nil "envs%i_l+~d" il1)
                  (format nil "envs%j_l+~a" r1) "j")
    (f90-g1e-emit fout op2 t2 s2 (format nil "envs%i_l+~d" il2)
                  (format nil "envs%j_l+~a" r2) "i")
    (format fout "      do ix = 0, envs%g_size*3 - 1~%")
    (format fout "         g(g~a+ix) = g(g~a+ix) + g(g~a+ix)~%" acc-a acc-a acc-b)
    (format fout "      end do~%")))

;;; The origin vector an RC/Ri/Rj/Rk operation is measured from.  The C
;;; emits these as a local double dr<symb>[3]; declaring and filling them is
;;; not optional -- an RC operation is a call to cint_x1?_1e taking that
;;; vector, and leaving it out is how the first version of this backend
;;; silently produced a plain overlap for int1e_r.
(defun dump-declare-dri-for-rc-f90 (fout i-ops symb)
  (flet ((emit-decl (from)
           (format fout "      real(dp) :: dr~a(0:2)~%" symb)
           (format fout "      dr~a = ~a~%" symb from)))
    (cond ((intersection '(rc xc yc zc) i-ops)
           (emit-decl (format nil "envs%r~a - envs%env(PTR_COMMON_ORIG:PTR_COMMON_ORIG+2)" symb)))
          ((intersection '(ri xi yi zi) i-ops)
           (emit-decl (format nil "envs%r~a - envs%ri" symb)))
          ((intersection '(rj xj yj zj) i-ops)
           (emit-decl (format nil "envs%r~a - envs%rj" symb)))
          ((intersection '(rk xk yk zk) i-ops)
           (emit-decl (format nil "envs%r~a - envs%rk" symb)))
          ((intersection '(rl xl yl zl) i-ops)
           (emit-decl (format nil "envs%r~a - envs%rl" symb))))))

;;; Fortran wants declarations before executable statements, so the two halves
;;; are emitted separately rather than interleaved the way the C does it.
(defun dri-decl-f90 (fout i-ops symb)
  (when (intersection '(rc xc yc zc ri xi yi zi rj xj yj zj rk xk yk zk rl xl yl zl) i-ops)
    (format fout "      real(dp) :: dr~a(0:2)~%" symb)))

(defun dri-init-f90 (fout i-ops symb)
  (cond ((intersection '(rc xc yc zc) i-ops)
         (format fout "      dr~a = envs%r~a - envs%env(PTR_COMMON_ORIG:PTR_COMMON_ORIG+2)~%" symb symb))
        ((intersection '(ri xi yi zi) i-ops)
         (format fout "      dr~a = envs%r~a - envs%ri~%" symb symb))
        ((intersection '(rj xj yj zj) i-ops)
         (format fout "      dr~a = envs%r~a - envs%rj~%" symb symb))
        ((intersection '(rk xk yk zk) i-ops)
         (format fout "      dr~a = envs%r~a - envs%rk~%" symb symb))
        ((intersection '(rl xl yl zl) i-ops)
         (format fout "      dr~a = envs%r~a - envs%rl~%" symb symb))))

;;; The GIAO constants.  A description carrying `g` contributes a factor of
;;; (r_i - r_j) -- or (r_k - r_l) -- per occurrence, and the C precomputes the
;;; 3^n products into a local c[] that the cell converter then indexes.  Two
;;; halves again, because Fortran wants declarations before executable
;;; statements.
(defun giao-decl-f90 (fout bra ket &optional (kbra '()) (kket '()))
  (let ((n-ij (count 'g (append bra ket)))
        (n-kl (count 'g (append kbra kket))))
    (when (> n-ij 0) (format fout "      real(dp) :: rirj(0:2)~%"))
    (when (> n-kl 0) (format fout "      real(dp) :: rkrl(0:2)~%"))
    (when (> (+ n-ij n-kl) 0)
      (format fout "      real(dp) :: c(0:~a)~%" (1- (expt 3 (+ n-ij n-kl)))))))

(defun giao-init-f90 (fout bra ket &optional (kbra '()) (kket '()))
  (let ((n-ij (count 'g (append bra ket)))
        (n-kl (count 'g (append kbra kket))))
    (when (> n-ij 0)
      (format fout "      rirj = envs%ri - envs%rj~%"))
    (when (> n-kl 0)
      (format fout "      rkrl = envs%rk - envs%rl~%"))
    (when (> (+ n-ij n-kl) 0)
      (loop
        for i upto (1- (expt 3 (+ n-ij n-kl))) do
        (format fout "      c(~a) = 1.0_dp" i)
        (loop
          for j from (+ n-ij n-kl -1) downto n-kl
          and res = i then (multiple-value-bind (int res) (floor res (expt 3 j))
                             (format fout " * rirj(~a)" int)
                             res))
        (loop
          for j from (1- n-kl) downto 0
          and res = (nth-value 1 (floor i (expt 3 n-kl)))
                    then (multiple-value-bind (int res) (floor res (expt 3 j))
                           (format fout " * rkrl(~a)" int)
                           res))
        (format fout "~%")))))

;;; The g-intermediate offsets.  A fourth-derivative integral has 32 of them,
;;; and one declaration line for those is 160 characters -- past Fortran's
;;; 132-column limit, which is an error rather than a warning.  Eight per line.
(defun emit-g-decls (fout ng)
  (let ((names (cons "g0" (loop for i in (range ng) collect (format nil "g~a" (1+ i))))))
    (loop for chunk on names by (lambda (l) (nthcdr 8 l)) do
          (format fout "      integer  :: ~{~a~^, ~}~%"
                  (subseq chunk 0 (min 8 (length chunk)))))))

;;; Wrapping the gout assembly.
;;;
;;; A fourth-derivative component can be eighty terms long, and Fortran's
;;; 132-column limit is an error rather than a warning.  Split at the term
;;; boundaries the cell converter puts in -- a space followed by + or -, which
;;; an exponent's minus never is, because that one follows an `e`.  The
;;; expression is unchanged and so is its left-to-right association; only the
;;; line breaks are new.
(defun split-gout-terms (str)
  (let ((res '()) (start 0) (n (length str)))
    (loop for i from 1 below n do
          (when (and (char= (char str i) #\Space)
                     (< (1+ i) n)
                     (member (char str (1+ i)) '(#\+ #\-)))
            (push (subseq str start i) res)
            (setf start i)))
    (push (subseq str start) res)
    (nreverse res)))

(defun emit-gout-line (fout lhs expr)
  (let ((line lhs)
        (has-term nil))
    (dolist (tm (split-gout-terms expr))
      (when (and has-term (> (+ (length line) (length tm)) 108))
        (format fout "~a &~%" line)
        (setf line "                       ")
        (setf has-term nil))
      (setf line (concatenate 'string line tm))
      (setf has-term t))
    (format fout "~a~%" line)))

;;; ---------------------------------------------------------------- gout

;;; The rinv family: nuc, rinv and nabla-rinv.  These have a Rys quadrature in
;;; them, so their gout sums over the roots and their recursion is the G2E_*
;;; one with k and l pinned at zero -- which is why D8 could not emit them and
;;; D9 can: the eight two-electron operations arrived with the relativistic
;;; work.
(defun gen-f90-gout1e-rinv (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore bra-k ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (ng (num-g-intermediates tot-bits op i-len j-len)))

      (setf *f90-gptr* (if (<= (1+ ng) *f90-gptr-max*) (1+ ng) nil))
      (format fout "   subroutine CINTgout1e_~a(gout, g, idx, envs, gout_empty)~%" intname)
      (format fout "      real(dp), intent(inout) :: gout(0:*)~%")
      (format fout "      real(dp), intent(inout), target :: g(0:)~%")
      (format fout "      integer,  intent(in)    :: idx(0:*)~%")
      (format fout "      type(cint_env_vars), intent(in) :: envs~%")
      (format fout "      integer,  intent(in)    :: gout_empty~%")
      (format fout "      integer  :: nf, nroots, ix, iy, iz, i, n~%")
      (emit-g-decls fout ng)
      (when *f90-gptr* (emit-gptr-decls fout *f90-gptr*))
      (format fout "      real(dp) :: s(0:~a)~%" (1- (expt 3 tot-bits)))
      (dri-decl-f90 fout bra-i "i")
      (dri-decl-f90 fout (append op ket-j) "j")
      (giao-decl-f90 fout bra-i (append op ket-j))
      (format fout "~%")
      (dri-init-f90 fout bra-i "i")
      (dri-init-f90 fout (append op ket-j) "j")
      (giao-init-f90 fout bra-i (append op ket-j))
      (format fout "      nf = envs%nf~%")
      (format fout "      nroots = envs%nrys_roots~%")
      (format fout "      g0 = 0~%")
      (loop for i in (range ng) do
            (format fout "      g~a = g~a + envs%g_size * 3~%" (1+ i) i))

      (dump-combo-braket fout (f90-fmt-i2-fn "envs%j_l" "0" "0")
                         (f90-fmt-op2-fn "envs%j_l" "0" "0")
                         (f90-fmt-j2-fn "envs%j_l" "0" "0")
                         i-rev op-rev j-rev 0)
      ;; deriv-max 0: always the accumulate loop, never the unrolled switch
      (when *f90-gptr* (emit-gptr-binds fout *f90-gptr*))
      (dump-s-2e-f90 fout tot-bits 0)
      (gen-f90-block-with-empty fout flat-script)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout1e_~a~%~%" intname)
      (setf *f90-gptr* nil))))

(defun gen-f90-gout1e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore bra-k ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (ng (num-g-intermediates tot-bits op i-len j-len)))

      (setf *f90-gptr* (if (<= (1+ ng) *f90-gptr-max*) (1+ ng) nil))
      (format fout "   subroutine CINTgout1e_~a(gout, g, idx, envs, gout_empty)~%" intname)
      (format fout "      real(dp), intent(inout) :: gout(0:*)~%")
      ;; g is INOUT: the G1E_* operations build intermediates inside it.
      (format fout "      real(dp), intent(inout), target :: g(0:)~%")
      (format fout "      integer,  intent(in)    :: idx(0:*)~%")
      (format fout "      type(cint_env_vars), intent(in) :: envs~%")
      (format fout "      integer,  intent(in)    :: gout_empty~%")
      (format fout "      integer  :: nf, ix, iy, iz, n~%")
      (emit-g-decls fout ng)
      (when *f90-gptr* (emit-gptr-decls fout *f90-gptr*))
      (format fout "      real(dp) :: s(0:~a)~%" (1- (expt 3 tot-bits)))
      (dri-decl-f90 fout bra-i "i")
      (dri-decl-f90 fout (append op ket-j) "j")
      (giao-decl-f90 fout bra-i (append op ket-j))
      (format fout "~%")
      (dri-init-f90 fout bra-i "i")
      (dri-init-f90 fout (append op ket-j) "j")
      (giao-init-f90 fout bra-i (append op ket-j))
      (format fout "      nf = envs%nf~%")
      (format fout "      g0 = 0~%")
      (loop for i in (range ng) do
            (format fout "      g~a = g~a + envs%g_size * 3~%" (1+ i) i))

      (dump-combo-braket fout (f90-fmt-i-fn) (f90-fmt-op-fn) (f90-fmt-j-fn)
                         i-rev op-rev j-rev 0)

      (when *f90-gptr* (emit-gptr-binds fout *f90-gptr*))
      (format fout "~%      do n = 0, nf - 1~%")
      (format fout "         ix = idx(0+n*3)~%")
      (format fout "         iy = idx(1+n*3)~%")
      (format fout "         iz = idx(2+n*3)~%")
      (dump-s-for-nroots-f90 fout tot-bits 1)
      (gen-f90-block-with-empty fout flat-script)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout1e_~a~%~%" intname)
      (setf *f90-gptr* nil))))

;;; -------------------------------------------------------------- drivers

;;; The optimizer entry point.  One per integral, as in the C, and each is a
;;; three-line wrapper around the builder for its arity -- which is all the
;;; C's 212 of them are too.  The ng it passes is the same literal the
;;; integral's own entry points pass, so an optimizer and the integral it
;;; belongs to always agree on what they were built for; the driver checks.
(defun gen-f90-optimizer (fout intname ngdef builder)
  (format fout "   subroutine ~a_optimizer(opt, atm, natm, bas, nbas, env)~%" intname)
  (format fout "      type(cint_opt_t), intent(inout) :: opt~%")
  (format fout "      integer,  intent(in) :: natm, nbas~%")
  (format fout "      integer,  intent(in) :: atm(0:), bas(0:)~%")
  (format fout "      real(dp), intent(in) :: env(0:)~%")
  (format fout "      integer, parameter :: ng(0:7) = ~a~%" ngdef)
  (format fout "      call ~a(opt, ng, atm, natm, bas, nbas, env)~%" builder)
  (format fout "   end subroutine ~a_optimizer~%~%" intname))

(defun gen-f90-int1e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore bra-k ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (raw-script (eval-int raw-infix))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (ts (car raw-script))
           (sf (cadr raw-script))
           (goutinc (length flat-script))
           (e1comps (if (eql sf 'sf) 1 4))
           (tensors (if (eql sf 'sf) goutinc (/ goutinc 4)))
           (rinv? (or (member 'nuc raw-infix) (member 'rinv raw-infix)
                      (member 'nabla-rinv raw-infix)))
           (int1e-type (cond ((member 'nuc raw-infix) 2)
                             ((or (member 'rinv raw-infix)
                                  (member 'nabla-rinv raw-infix)) 1)
                             (t 0)))
           (ngdef (if rinv?
                      (format nil "[~d, ~d, 0, 0, ~d, ~d, 0, ~d]"
                              (if (intersection *act-left-right* op) (1+ i-len) i-len)
                              (+ op-len j-len) tot-bits e1comps tensors)
                      (format nil "[~d, ~d, 0, 0, ~d, ~d, 1, ~d]"
                              i-len (+ op-len j-len) tot-bits e1comps tensors)))
           (fac (factor-of raw-infix)))

      (format fout "   ! <~{~a ~}i|~{~a ~}|~{~a ~}j>~%" bra-i op ket-j)
      (gen-f90-optimizer fout intname ngdef "cint_all_1e_optimizer")
      (if rinv?
          (gen-f90-gout1e-rinv fout intname raw-infix flat-script)
          (gen-f90-gout1e fout intname raw-infix flat-script))

      (dolist (form '(("cart" . "C2S_CART_1E") ("sph" . "C2S_SPH_1E")))
        (format fout "   function ~a_~a(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%"
                intname (car form))
        (format fout "         result(has_value)~%")
        (format fout "      real(dp), intent(inout) :: out(0:)~%")
        (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), atm(0:), bas(0:)~%")
        (format fout "      integer,  intent(in)    :: natm, nbas~%")
        (format fout "      real(dp), intent(in)    :: env(0:)~%")
        (format fout "      type(cint_ws), intent(inout) :: ws~%")
        (format fout "      logical :: has_value~%")
        (format fout "      type(cint_env_vars) :: envs~%")
        (format fout "      integer, parameter :: ng(0:7) = ~a~%" ngdef)
        (format fout "~%")
        (format fout "      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)~%")
        (format fout "      envs%f_gout => CINTgout1e_~a~%" intname)
        (unless (= fac 1)
          (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
        (format fout "      has_value = cint_1e_drv(out, dims, envs, ws, ~a, ~d)~%"
                (cdr form) int1e-type)
        (format fout "   end function ~a_~a~%~%" intname (car form)))

      ;; _spinor.  Same envs, same gout, same primitive loop -- only the
      ;; transform on the way out and the output type differ, which is the
      ;; whole reason the spinor forms are nearly free once the runtime is
      ;; there.  (ts, sf) pick which of the four c2s_s[fi]_1e[i] runs, as the
      ;; two booleans cint_1e_spinor_drv takes.
      (format fout "   function ~a_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%" intname)
      (format fout "         result(has_value)~%")
      (format fout "      complex(dp), intent(inout) :: out(0:)~%")
      (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas~%")
      (format fout "      integer,  target        :: atm(0:), bas(0:)~%")
      (format fout "      real(dp), target        :: env(0:)~%")
      (format fout "      type(cint_ws), intent(inout) :: ws~%")
      (format fout "      logical :: has_value~%")
      (format fout "      type(cint_env_vars) :: envs~%")
      (format fout "      integer, parameter :: ng(0:7) = ~a~%~%" ngdef)
      (format fout "      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)~%")
      (format fout "      envs%f_gout => CINTgout1e_~a~%" intname)
      (unless (= fac 1)
        (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
      (format fout "      has_value = cint_1e_spinor_drv(out, dims, envs, ws, ~a, ~d)~%"
              (f90-c2s-flags sf ts) int1e-type)
      (format fout "   end function ~a_spinor~%~%" intname))))

;;; ---------------------------------------------------------------- driver


;;; ================================================================ 2e
;;;
;;; Everything below is the four-centre counterpart of the one-centre code
;;; above.  The symbolic walk is again gen-code.cl's, unchanged.

;;; The s-accumulation.  The C switches on nrys_roots and unrolls one to four
;;; roots into a single expression, falling back to a zero-then-accumulate
;;; loop above that.  Both forms are reproduced rather than collapsed into the
;;; loop: they agree in IEEE only because adding a term to a fresh zero is
;;; exact, and that is an argument worth having made once, in the open, rather
;;; than assumed at every call site.
(defun dump-s-2e-f90 (fout tot-bits &optional (deriv-max 2))
  (format fout "~%      do n = 0, nf - 1~%")
  (format fout "         ix = idx(0+n*3)~%")
  (format fout "         iy = idx(1+n*3)~%")
  (format fout "         iz = idx(2+n*3)~%")
  (if (< tot-bits deriv-max)
    (progn
      (format fout "         select case (nroots)~%")
      (loop for i from 1 to 4 do
            (format fout "         case (~a)~%" i)
            (dump-s-for-nroots-f90 fout tot-bits i "            "))
      (format fout "         case default~%")
      (format fout "            s = 0.0_dp~%")
      (format fout "            do i = 0, nroots - 1~%")
      (dump-s-loop-f90 fout tot-bits)
      (format fout "            end do~%")
      (format fout "         end select~%"))
    (progn
      (format fout "         s = 0.0_dp~%")
      (format fout "         do i = 0, nroots - 1~%")
      (dump-s-loop-f90 fout tot-bits)
      (format fout "         end do~%"))))

(defun dump-s-loop-f90 (fout tot-bits)
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "               s(~a) = s(~a) + ~a * ~a * ~a~%"
              i i (gref xbin "ix+i") (gref ybin "iy+i") (gref zbin "iz+i")))))

;;; The G2E_* split, exactly parallel to the G1E_* one: calls, pointer shifts,
;;; and which origin vector an x-operation measures from.
(defparameter *f90-2e-call*
  '(("D_" . "cint_nabla1~a_2e") ("R0" . "cint_x1~a_2e") ("RC" . "cint_x1~a_2e")))
(defparameter *f90-2e-shift*
  '(("R_" . "g_stride_~a")))
(defparameter *f90-2e-origin*
  '(("R0" . "envs%r~a") ("RC" . "dr~a")))

(defun f90-g2e-emit (fout op tgt src il jl kl ll centre)
  (let ((call  (cdr (assoc op *f90-2e-call*  :test #'string=)))
        (shift (cdr (assoc op *f90-2e-shift* :test #'string=)))
        (orig  (cdr (assoc op *f90-2e-origin* :test #'string=)))
        (c (string-downcase centre)))
    (cond
      (call
        (if orig
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src (format nil orig c) il jl kl ll)
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src il jl kl ll)))
      (shift
        (format fout "      g~a = g~a + envs%~a~%" tgt src (format nil shift c)))
      (t (error "unmapped G2E op ~a" op)))))

;;; Closures for the k/l walk and the i/j walk.  Which angular momenta are
;;; held fixed, and at what, is an argument rather than baked into a format
;;; string -- that is the whole difference between the four-, three- and
;;; two-centre cases, and the C spells it by writing literal zeros into three
;;; near-identical copies.
(defun f90-fmt-k-fn (ifix jfix lfix)
  (lambda (fout op tgt src lshift)
    (f90-g2e-emit fout op tgt src ifix jfix
                  (format nil "envs%k_l+~a" lshift) lfix "k")))

(defun f90-fmt-l-fn (ifix jfix)
  (lambda (fout op tgt src klen right)
    (f90-g2e-emit fout op tgt src ifix jfix
                  (format nil "envs%k_l+~a" klen)
                  (format nil "envs%l_l+~a" right) "l")))

(defun f90-fmt-i2-fn (jbase kfix lfix)
  (lambda (fout op tgt src lshift)
    (f90-g2e-emit fout op tgt src
                  (format nil "envs%i_l+~a" lshift) jbase kfix lfix "i")))

(defun f90-fmt-j2-fn (jbase kfix lfix)
  (lambda (fout op tgt src ilen right)
    (f90-g2e-emit fout op tgt src
                  (format nil "envs%i_l+~d" ilen)
                  (format nil "~a+~a" jbase right) kfix lfix "j")))

(defun f90-fmt-op2-fn (jbase kfix lfix)
  (lambda (fout op1 t1 s1 il1 r1 op2 t2 s2 il2 r2 acc-a acc-b)
    (f90-g2e-emit fout op1 t1 s1 (format nil "envs%i_l+~d" il1)
                  (format nil "~a+~a" jbase r1) kfix lfix "j")
    (f90-g2e-emit fout op2 t2 s2 (format nil "envs%i_l+~a" il2)
                  (format nil "~a+~a" jbase r2) kfix lfix "i")
    (format fout "      do ix = 0, envs%g_size*3 - 1~%")
    (format fout "         g(g~a+ix) = g(g~a+ix) + g(g~a+ix)~%" acc-a acc-a acc-b)
    (format fout "      end do~%")))

;;; The declaration block every gout2e-family kernel opens with.  Shared
;;; because the three arities differ only in which dr vectors they declare.
(defun gen-f90-gout2e-head (fout intname ng tot-bits bra-i op ket-j bra-k ket-l)
  (format fout "   subroutine CINTgout2e_~a(gout, g, idx, envs, gout_empty)~%" intname)
  (format fout "      real(dp), intent(inout) :: gout(0:*)~%")
  (format fout "      real(dp), intent(inout), target :: g(0:)~%")
  (format fout "      integer,  intent(in)    :: idx(0:*)~%")
  (format fout "      type(cint_env_vars), intent(in) :: envs~%")
  (format fout "      integer,  intent(in)    :: gout_empty~%")
  (format fout "      integer  :: nf, nroots, ix, iy, iz, i, n~%")
  (emit-g-decls fout ng)
  (when *f90-gptr* (emit-gptr-decls fout *f90-gptr*))
  (format fout "      real(dp) :: s(0:~a)~%" (1- (expt 3 tot-bits)))
  (dri-decl-f90 fout bra-i "i")
  (dri-decl-f90 fout (append op ket-j) "j")
  (dri-decl-f90 fout bra-k "k")
  (dri-decl-f90 fout ket-l "l")
  (giao-decl-f90 fout bra-i (append op ket-j) bra-k ket-l)
  (format fout "~%")
  (dri-init-f90 fout bra-i "i")
  (dri-init-f90 fout (append op ket-j) "j")
  (dri-init-f90 fout bra-k "k")
  (dri-init-f90 fout ket-l "l")
  (giao-init-f90 fout bra-i (append op ket-j) bra-k ket-l)
  (format fout "      nf = envs%nf~%")
  (format fout "      nroots = envs%nrys_roots~%")
  (format fout "      g0 = 0~%")
  (loop for i in (range ng) do
        (format fout "      g~a = g~a + envs%g_size * 3~%" (1+ i) i)))

;;; One entry point.  `body` is the driver call, `initfn` the envs setup, and
;;; `spinor` decides the output type -- everything else is the same six
;;; declarations for all three arities and both angular forms.
(defun gen-f90-2e-entry (fout intname suffix ngdef fac spinor body initfn
                         &optional (post nil))
  (format fout "   function ~a_~a(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%"
          intname suffix)
  (format fout "         result(has_value)~%")
  (if spinor
      (format fout "      complex(dp), intent(inout) :: out(0:)~%")
      (format fout "      real(dp), intent(inout) :: out(0:)~%"))
  (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas~%")
  (format fout "      integer,  target        :: atm(0:), bas(0:)~%")
  (format fout "      real(dp), target        :: env(0:)~%")
  (format fout "      type(cint_ws), intent(inout) :: ws~%")
  (format fout "      logical :: has_value~%")
  (format fout "      type(cint_env_vars) :: envs~%")
  (format fout "      integer, parameter :: ng(0:7) = ~a~%~%" ngdef)
  (format fout "      call ~a(envs, ng, shls, atm, natm, bas, nbas, env)~%" initfn)
  (format fout "      envs%f_gout => CINTgout2e_~a~%" intname)
  (unless (= fac 1)
    (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
  (when post (format fout "~a" post))
  (format fout "      has_value = ~a~%" body)
  (format fout "   end function ~a_~a~%~%" intname suffix))

(defun gen-f90-gout2e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (l-rev (reverse (effect-keys ket-l)))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (l-len (length l-rev))
           (tot-bits (+ i-len j-len op-len k-len l-len))
           (ng (num-g-intermediates tot-bits op i-len j-len))
           (act-lr (intersection *act-left-right* op)))

      (setf *f90-gptr* (if (<= (1+ ng) *f90-gptr-max*) (1+ ng) nil))
      (gen-f90-gout2e-head fout intname ng tot-bits bra-i op ket-j bra-k ket-l)

      ;; the k/l walk first, then i/j, as the C does -- the order matters
      ;; because the second walk reads intermediates the first one wrote
      (let ((ifix (if act-lr (format nil "envs%i_l+~a" (1+ i-len))
                             (format nil "envs%i_l+~a" i-len)))
            (jfix (if act-lr (format nil "envs%j_l+~a" (1+ j-len))
                             (format nil "envs%j_l+~a" (+ op-len j-len)))))
        (dump-combo-braket fout (f90-fmt-k-fn ifix jfix "envs%l_l")
                           (lambda (&rest a) (declare (ignore a))
                             (error "no operator on the k/l side"))
                           (f90-fmt-l-fn ifix jfix)
                           k-rev '() l-rev 0))
      (dump-combo-braket fout (f90-fmt-i2-fn "envs%j_l" "envs%k_l" "envs%l_l")
                         (f90-fmt-op2-fn "envs%j_l" "envs%k_l" "envs%l_l")
                         (f90-fmt-j2-fn "envs%j_l" "envs%k_l" "envs%l_l")
                         i-rev op-rev j-rev (+ k-len l-len))

      (when *f90-gptr* (emit-gptr-binds fout *f90-gptr*))
      (dump-s-2e-f90 fout tot-bits)
      (gen-f90-block-with-empty fout flat-script)
      (setf *f90-gptr* nil)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout2e_~a~%~%" intname))))

;;; Which of the eight complex c2s transforms a spinor stage uses, as the two
;;; booleans cint_2e_spinor_drv takes: spin-included, and carrying a factor i.
(defun f90-c2s-flags (sf ts)
  (format nil "~a, ~a"
          (if (eql sf 'sf) ".false." ".true.")
          (if (eql ts 'ts) ".false." ".true.")))

(defun gen-f90-int4c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (l-rev (reverse (effect-keys ket-l)))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (l-len (length l-rev))
           (tot-bits (+ i-len j-len op-len k-len l-len))
           (raw-script (eval-int raw-infix))
           (ts1 (car raw-script))
           (sf1 (cadr raw-script))
           (ts2 (caddr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (if (intersection *act-left-right* op)
                      (format nil "[~d, ~d, ~d, ~d, ~d, ~d, ~d, ~d]"
                              (1+ i-len) (1+ j-len) k-len l-len tot-bits e1comps e2comps tensors)
                      (format nil "[~d, ~d, ~d, ~d, ~d, ~d, ~d, ~d]"
                              i-len (+ op-len j-len) k-len l-len tot-bits e1comps e2comps tensors)))
           (fac (factor-of raw-infix))
           ;; the C flips the sign when both electrons are spin-included and
           ;; time-antisymmetric; two factors of i on the same integral
           (flip (and (eql sf1 'si) (eql ts1 'tas) (eql sf2 'si) (eql ts2 'tas))))

      (format fout "   ! (~{~a ~}i ~{~a ~}j|~{~a ~}|~{~a ~}k ~{~a ~}l)~%"
              bra-i ket-j op bra-k ket-l)
      (gen-f90-optimizer fout intname ngdef "cint_all_2e_optimizer")
      (gen-f90-gout2e fout intname raw-infix flat-script)

      ;; The C flips the sign in _cart and _sph and only there: two factors
      ;; of i on the same integral, which the spinor transform carries itself.
      (let ((flipline (when flip
                        (format nil "      ! two factors of i on the same integral~%      envs%common_factor = -envs%common_factor~%"))))
        (dolist (form '(("cart" . "C2S_CART_2E1") ("sph" . "C2S_SPH_2E1")))
          (gen-f90-2e-entry fout intname (car form) ngdef fac nil
                            (format nil "cint_2e_drv(out, dims, envs, ws, ~a)" (cdr form))
                            "cint_init_int2e_envvars" flipline)))
      (gen-f90-2e-entry fout intname "spinor" ngdef fac t
                        (format nil "cint_2e_spinor_drv(out, dims, envs, ws, ~a, ~a)"
                                (f90-c2s-flags sf1 ts1) (f90-c2s-flags sf2 ts2))
                        "cint_init_int2e_envvars"))))


;;; ============================================================ 3c2e, 2c2e
;;;
;;; Both are the four-centre walk with indices held at zero: 3c2e drops l,
;;; 2c2e drops j and l as well.  The C says that by writing literal zeros
;;; into the format strings; here the fixed text is an argument, so the same
;;; closures serve all three arities.

(defun gen-f90-gout3c2e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (ng (num-g-intermediates tot-bits op i-len j-len))
           (act-lr (intersection *act-left-right* op)))

      (gen-f90-gout2e-head fout intname ng tot-bits bra-i op ket-j bra-k '())

      ;; the k walk, with l pinned at zero
      (let ((ifix (if act-lr (format nil "envs%i_l+~a" (1+ i-len))
                             (format nil "envs%i_l+~a" i-len)))
            (jfix (if act-lr (format nil "envs%j_l+~a" (1+ j-len))
                             (format nil "envs%j_l+~a" (+ op-len j-len)))))
        (dump-combo-braket fout (f90-fmt-k-fn ifix jfix "0")
                           (lambda (&rest a) (declare (ignore a))
                             (error "no operator on the k side"))
                           (lambda (&rest a) (declare (ignore a))
                             (error "no l index in a three-centre integral"))
                           k-rev '() '() 0))
      (dump-combo-braket fout (f90-fmt-i2-fn "envs%j_l" "envs%k_l" "0")
                         (f90-fmt-op2-fn "envs%j_l" "envs%k_l" "0")
                         (f90-fmt-j2-fn "envs%j_l" "envs%k_l" "0")
                         i-rev op-rev j-rev k-len)

      (dump-s-2e-f90 fout tot-bits)
      (gen-f90-block-with-empty fout flat-script)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout2e_~a~%~%" intname))))

(defun gen-f90-gout2c2e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-j ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len op-len k-len))
           (ng (num-g-intermediates tot-bits op i-len 0))
           (act-lr (intersection *act-left-right* op)))

      (gen-f90-gout2e-head fout intname ng tot-bits bra-i op '() bra-k '())

      (let ((ifix (if act-lr (format nil "envs%i_l+~a" (1+ i-len))
                             (format nil "envs%i_l+~a" i-len))))
        (dump-combo-braket fout (f90-fmt-k-fn ifix "0" "0")
                           (lambda (&rest a) (declare (ignore a))
                             (error "no operator on the k side"))
                           (lambda (&rest a) (declare (ignore a))
                             (error "no l index in a two-centre integral"))
                           k-rev '() '() 0))
      (dump-combo-braket fout (f90-fmt-i2-fn "0" "envs%k_l" "0")
                         (f90-fmt-op2-fn "0" "envs%k_l" "0")
                         (f90-fmt-j2-fn "0" "envs%k_l" "0")
                         i-rev op-rev '() k-len)

      (dump-s-2e-f90 fout tot-bits)
      (gen-f90-block-with-empty fout flat-script)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout2e_~a~%~%" intname))))

(defun gen-f90-int3c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (ts1 (car raw-script))
           (sf1 (cadr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (if (intersection *act-left-right* op)
                      (format nil "[~d, ~d, ~d, 0, ~d, ~d, ~d, ~d]"
                              (1+ i-len) (1+ j-len) k-len tot-bits e1comps e2comps tensors)
                      (format nil "[~d, ~d, ~d, 0, ~d, ~d, ~d, ~d]"
                              i-len (+ op-len j-len) k-len tot-bits e1comps e2comps tensors)))
           (fac (factor-of raw-infix)))

      (format fout "   ! (~{~a ~}i ~{~a ~}j|~{~a ~}|~{~a ~}k)~%" bra-i ket-j op bra-k)
      (gen-f90-optimizer fout intname ngdef "cint_all_3c2e_optimizer")
      (gen-f90-gout3c2e fout intname raw-infix flat-script)
      (dolist (form '(("cart" . "C2S_CART_3C2E1") ("sph" . "C2S_SPH_3C2E1")))
        (gen-f90-2e-entry fout intname (car form) ngdef fac nil
                          (format nil "cint_3c2e_drv(out, dims, envs, ws, ~a)" (cdr form))
                          "cint_init_int3c2e_envvars"))
      (gen-f90-2e-entry fout intname "spinor" ngdef fac t
                        (format nil "cint_3c2e_spinor_drv(out, dims, envs, ws, ~a, .false.)"
                                (f90-c2s-flags sf1 ts1))
                        "cint_init_int3c2e_envvars"))))

(defun gen-f90-int2c2e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-j ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (sf1 (cadr raw-script))
           (sf2 (cadddr raw-script))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (goutinc (length flat-script))
           (e1comps (if (eql sf1 'sf) 1 4))
           (e2comps (if (eql sf2 'sf) 1 4))
           (tensors (cond ((and (eql sf1 'sf) (eql sf2 'sf)) goutinc)
                          ((and (eql sf1 'si) (eql sf2 'si)) (/ goutinc 16))
                          (t (/ goutinc 4))))
           (ngdef (if (intersection *act-left-right* op)
                      (format nil "[~d, 0, ~d, 0, ~d, ~d, ~d, ~d]"
                              (1+ i-len) k-len tot-bits e1comps e2comps tensors)
                      (format nil "[~d, ~d, ~d, 0, ~d, ~d, ~d, ~d]"
                              i-len op-len k-len tot-bits e1comps e2comps tensors)))
           (fac (factor-of raw-infix)))

      (format fout "   ! (~{~a ~}i|~{~a ~}|~{~a ~}k)~%" bra-i op bra-k)
      (gen-f90-optimizer fout intname ngdef "cint_all_2c2e_optimizer")
      (gen-f90-gout2c2e fout intname raw-infix flat-script)
      (dolist (form '(("cart" . "C2S_CART_2C2E1") ("sph" . "C2S_SPH_2C2E1")))
        (gen-f90-2e-entry fout intname (car form) ngdef fac nil
                          (format nil "cint_2c2e_drv(out, dims, envs, ws, ~a)" (cdr form))
                          "cint_init_int2c2e_envvars"))
      ;; No _spinor: CINT2c2e_spinor_drv is an unimplemented stub upstream
      ;; (cint2c2e.c:294-298), and the C's own generator emits a function that
      ;; prints "not implemented" and returns.  Inherited rather than invented.
      )))



;;; ============================================================== 1e grids
;;;
;;; The grid index is innermost and blocked, so both the s-accumulation and
;;; the gout assembly carry an extra loop over the block -- which is what the
;;; `with-grids` flag in the cell converter is for, and why these get their
;;; own emitters rather than a flag on the existing ones.

(defun f90-cell-converter-grids (cell fout &optional (with-grids nil))
  (declare (ignore with-grids))
  (f90-cell-converter cell fout t))

(defun gen-f90-block-grids (fout flat-script accumulate)
  (let ((assemb (to-c-code-string fout #'f90-cell-converter-grids flat-script t))
        (comp (length flat-script)))
    (loop for s in assemb
          for gid from 0 do
          (if accumulate
              (emit-gout-line fout
                              (format nil "               gout(ig+bgrids*(n*~a+~a)) = gout(ig+bgrids*(n*~a+~a)) +"
                                      comp gid comp gid)
                              (format nil " (~a )" s))
              (emit-gout-line fout
                              (format nil "               gout(ig+bgrids*(n*~a+~a)) =" comp gid) s)))))

(defparameter *f90-grids-call*
  '(("D_" . "cint_nabla1~a_grids") ("R0" . "cint_x1~a_grids") ("RC" . "cint_x1~a_grids")))
(defparameter *f90-grids-shift*
  '(("R_" . "g_stride_~a")))

(defun f90-g1e-grids-emit (fout op tgt src il jl centre)
  (let ((call  (cdr (assoc op *f90-grids-call*  :test #'string=)))
        (shift (cdr (assoc op *f90-grids-shift* :test #'string=)))
        (orig  (cdr (assoc op *f90-1e-origin*   :test #'string=)))
        (c (string-downcase centre)))
    (cond
      (call
        (if orig
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src (format nil orig c) il jl)
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src il jl)))
      (shift
        (format fout "      g~a = g~a + envs%~a~%" tgt src (format nil shift c)))
      (t (error "unmapped G1E_GRIDS op ~a" op)))))

(defun f90-fmt-gi-fn ()
  (lambda (fout op tgt src lshift)
    (f90-g1e-grids-emit fout op tgt src
                        (format nil "envs%i_l+~a" lshift) "envs%j_l" "i")))

(defun f90-fmt-gj-fn ()
  (lambda (fout op tgt src ilen right)
    (f90-g1e-grids-emit fout op tgt src
                        (format nil "envs%i_l+~d" ilen)
                        (format nil "envs%j_l+~a" right) "j")))

;;; The mixed one: the C writes G1E_GRIDS_?J for the first and plain G1E_?I
;;; for the second, which is not a typo -- the operator acting between bra and
;;; ket is applied with the ordinary 1e recursion.
(defun f90-fmt-gop-fn ()
  (lambda (fout op1 t1 s1 il1 r1 op2 t2 s2 il2 r2 acc-a acc-b)
    (f90-g1e-grids-emit fout op1 t1 s1 (format nil "envs%i_l+~d" il1)
                        (format nil "envs%j_l+~a" r1) "j")
    (f90-g1e-emit fout op2 t2 s2 (format nil "envs%i_l+~d" il2)
                  (format nil "envs%j_l+~a" r2) "i")
    (format fout "      do ix = 0, envs%g_size*3 - 1~%")
    (format fout "         g(g~a+ix) = g(g~a+ix) + g(g~a+ix)~%" acc-a acc-a acc-b)
    (format fout "      end do~%")))

(defun gen-f90-gout1e-grids (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore bra-k ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (comp (expt 3 tot-bits))
           (ng (num-g-intermediates tot-bits op i-len j-len)))

      (setf *f90-gptr* (if (<= (1+ ng) *f90-gptr-max*) (1+ ng) nil))
      (format fout "   subroutine CINTgout1e_~a(gout, g, idx, envs, gout_empty)~%" intname)
      (format fout "      real(dp), intent(inout) :: gout(0:*)~%")
      (format fout "      real(dp), intent(inout), target :: g(0:)~%")
      (format fout "      integer,  intent(in)    :: idx(0:*)~%")
      (format fout "      type(cint_env_vars), intent(in) :: envs~%")
      (format fout "      integer,  intent(in)    :: gout_empty~%")
      (format fout "      integer  :: ngrids, bgrids, nroots, nf, ix, iy, iz, n, i, ig~%")
      (emit-g-decls fout ng)
      (when *f90-gptr* (emit-gptr-decls fout *f90-gptr*))
      (format fout "      real(dp) :: s(0:GRID_BLKSIZE*~a - 1)~%" comp)
      (dri-decl-f90 fout bra-i "i")
      (dri-decl-f90 fout (append op ket-j) "j")
      (giao-decl-f90 fout bra-i (append op ket-j))
      (format fout "~%")
      (dri-init-f90 fout bra-i "i")
      (dri-init-f90 fout (append op ket-j) "j")
      (giao-init-f90 fout bra-i (append op ket-j))
      (format fout "      ngrids = envs%nfl~%")
      (format fout "      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)~%")
      (format fout "      nroots = envs%nrys_roots~%")
      (format fout "      nf = envs%nf~%")
      (format fout "      g0 = 0~%")
      (loop for i in (range ng) do
            (format fout "      g~a = g~a + envs%g_size * 3~%" (1+ i) i))

      (dump-combo-braket fout (f90-fmt-gi-fn) (f90-fmt-gop-fn) (f90-fmt-gj-fn)
                         i-rev op-rev j-rev 0)

      (when *f90-gptr* (emit-gptr-binds fout *f90-gptr*))
      (format fout "~%      do n = 0, nf - 1~%")
      (format fout "         ix = idx(0+n*3)~%")
      (format fout "         iy = idx(1+n*3)~%")
      (format fout "         iz = idx(2+n*3)~%")
      (format fout "         do i = 0, ~a - 1~%" comp)
      (format fout "            do ig = 0, bgrids - 1~%")
      (format fout "               s(ig+i*GRID_BLKSIZE) = 0.0_dp~%")
      (format fout "            end do~%")
      (format fout "         end do~%")
      (format fout "         do i = 0, nroots - 1~%")
      (format fout "            do ig = 0, bgrids - 1~%")
      (dump-s-loop-grids-f90 fout tot-bits)
      (format fout "            end do~%")
      (format fout "         end do~%")
      (format fout "         if (gout_empty /= 0) then~%")
      (format fout "            do ig = 0, bgrids - 1~%")
      (gen-f90-block-grids fout flat-script nil)
      (format fout "            end do~%")
      (format fout "         else~%")
      (format fout "            do ig = 0, bgrids - 1~%")
      (gen-f90-block-grids fout flat-script t)
      (format fout "            end do~%")
      (format fout "         end if~%")
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout1e_~a~%~%" intname)
      (setf *f90-gptr* nil))))

(defun dump-s-loop-grids-f90 (fout tot-bits)
  (loop
    for i upto (1- (expt 3 tot-bits)) do
    (let* ((ybin (dec-to-ybin i))
           (zbin (dec-to-zbin i))
           (xbin (- (ash 1 tot-bits) 1 ybin zbin)))
      (format fout "               s(ig+GRID_BLKSIZE*~a) = s(ig+GRID_BLKSIZE*~a) &~%" i i)
      (format fout "                  + ~a * ~a &~%"
              (gref xbin "ix+ig+i*GRID_BLKSIZE") (gref ybin "iy+ig+i*GRID_BLKSIZE"))
      (format fout "                  * ~a~%" (gref zbin "iz+ig+i*GRID_BLKSIZE")))))

(defun gen-f90-int1e-grids (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore bra-k ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (op-rev (reverse (effect-keys op)))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (op-len (length op-rev))
           (tot-bits (+ i-len j-len op-len))
           (raw-script (eval-int raw-infix))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (ts (car raw-script))
           (sf (cadr raw-script))
           (goutinc (length flat-script))
           (e1comps (if (eql sf 'sf) 1 4))
           (tensors (if (eql sf 'sf) goutinc (/ goutinc 4)))
           (ngdef (if (intersection *act-left-right* op)
                      (format nil "[~d, ~d, 0, 0, ~d, ~d, 0, ~d]"
                              (1+ i-len) (+ op-len j-len) tot-bits e1comps tensors)
                      (format nil "[~d, ~d, 0, 0, ~d, ~d, 0, ~d]"
                              i-len (+ op-len j-len) tot-bits e1comps tensors)))
           (fac (factor-of raw-infix)))

      (format fout "   ! <~{~a ~}i|~{~a ~}grids|~{~a ~}j>~%" bra-i op ket-j)
      (gen-f90-optimizer fout intname ngdef "cint_all_1e_grids_optimizer")
      (gen-f90-gout1e-grids fout intname raw-infix flat-script)
      (dolist (form '(("cart" . "C2S_CART_1E") ("sph" . "C2S_SPH_1E")))
        (format fout "   function ~a_~a(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%"
                intname (car form))
        (format fout "         result(has_value)~%")
        (format fout "      real(dp), intent(inout) :: out(0:)~%")
        (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas~%")
        (format fout "      integer,  target        :: atm(0:), bas(0:)~%")
        (format fout "      real(dp), target        :: env(0:)~%")
        (format fout "      type(cint_ws), intent(inout) :: ws~%")
        (format fout "      logical :: has_value~%")
        (format fout "      type(cint_env_vars) :: envs~%")
        (format fout "      integer, parameter :: ng(0:7) = ~a~%~%" ngdef)
        (format fout "      call cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)~%")
        (format fout "      envs%f_gout => CINTgout1e_~a~%" intname)
        (unless (= fac 1)
          (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
        (format fout "      has_value = cint_1e_grids_drv(out, dims, envs, ws, ~a)~%" (cdr form))
        (format fout "   end function ~a_~a~%~%" intname (car form)))

      ;; _spinor.  As for the plain one-electron forms, the primitive loop and
      ;; the gout are shared and only the transform and the output type change.
      ;; The C reaches for c2s_sf_1e_grids or c2s_si_1e_grids and never for
      ;; either of the two _gridsi variants -- no grids integral it generates
      ;; carries an odd number of momentum factors -- but the driver takes the
      ;; same two booleans the one-electron one does, so (sf, ts) selects here
      ;; exactly as it does there.
      (format fout "   function ~a_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%" intname)
      (format fout "         result(has_value)~%")
      (format fout "      complex(dp), intent(inout) :: out(0:)~%")
      (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas~%")
      (format fout "      integer,  target        :: atm(0:), bas(0:)~%")
      (format fout "      real(dp), target        :: env(0:)~%")
      (format fout "      type(cint_ws), intent(inout) :: ws~%")
      (format fout "      logical :: has_value~%")
      (format fout "      type(cint_env_vars) :: envs~%")
      (format fout "      integer, parameter :: ng(0:7) = ~a~%~%" ngdef)
      (format fout "      call cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)~%")
      (format fout "      envs%f_gout => CINTgout1e_~a~%" intname)
      (unless (= fac 1)
        (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
      (format fout "      has_value = cint_1e_grids_spinor_drv(out, dims, envs, ws, ~a)~%"
              (f90-c2s-flags sf ts))
      (format fout "   end function ~a_spinor~%~%" intname))))


;;; ================================================================= 3c1e
;;;
;;; The three-centre one-electron family.  Its recursion is the ordinary 1e
;;; one -- G1E_?I/J/K, which the port already has, k variants included -- so
;;; the only new thing here is the k walk and a third fixed angular momentum
;;; in each call.

(defun f90-fmt-1e-k-fn (ifix jfix)
  (lambda (fout op tgt src lshift)
    (f90-g1e-emit-3 fout op tgt src ifix jfix
                    (format nil "envs%k_l+~a" lshift) "k")))

(defun f90-fmt-1e-i3-fn ()
  (lambda (fout op tgt src lshift)
    (f90-g1e-emit-3 fout op tgt src
                    (format nil "envs%i_l+~a" lshift) "envs%j_l" "envs%k_l" "i")))

(defun f90-fmt-1e-j3-fn ()
  (lambda (fout op tgt src ilen right)
    (f90-g1e-emit-3 fout op tgt src
                    (format nil "envs%i_l+~d" ilen)
                    (format nil "envs%j_l+~a" right) "envs%k_l" "j")))

(defun f90-fmt-1e-op3-fn ()
  (lambda (fout op1 t1 s1 il1 r1 op2 t2 s2 il2 r2 acc-a acc-b)
    (f90-g1e-emit-3 fout op1 t1 s1 (format nil "envs%i_l+~d" il1)
                    (format nil "envs%j_l+~a" r1) "envs%k_l" "j")
    (f90-g1e-emit-3 fout op2 t2 s2 (format nil "envs%i_l+~a" il2)
                    (format nil "envs%j_l+~a" r2) "envs%k_l" "i")
    (format fout "      do ix = 0, envs%g_size*3 - 1~%")
    (format fout "         g(g~a+ix) = g(g~a+ix) + g(g~a+ix)~%" acc-a acc-a acc-b)
    (format fout "      end do~%")))

;;; The same split as f90-g1e-emit, with the third angular momentum passed
;;; through rather than pinned.
(defun f90-g1e-emit-3 (fout op tgt src il jl kl centre)
  (let ((call  (cdr (assoc op *f90-1e-call*  :test #'string=)))
        (shift (cdr (assoc op *f90-1e-shift* :test #'string=)))
        (orig  (cdr (assoc op *f90-1e-origin* :test #'string=)))
        (c (string-downcase centre)))
    (cond
      (call
        (if orig
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src (format nil orig c) il jl kl)
            (format fout "      call ~a(g, g~a, g~a, ~a, ~a, ~a, envs)~%"
                    (format nil call c) tgt src il jl kl)))
      (shift
        (format fout "      g~a = g~a + envs%~a~%" tgt src (format nil shift c)))
      (t (error "unmapped G1E op ~a" op)))))

(defun gen-f90-gout3c1e (fout intname raw-infix flat-script)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (ng (num-g-intermediates tot-bits op i-len j-len)))

      (setf *f90-gptr* (if (<= (1+ ng) *f90-gptr-max*) (1+ ng) nil))
      (format fout "   subroutine CINTgout1e_~a(gout, g, idx, envs, gout_empty)~%" intname)
      (format fout "      real(dp), intent(inout) :: gout(0:*)~%")
      (format fout "      real(dp), intent(inout), target :: g(0:)~%")
      (format fout "      integer,  intent(in)    :: idx(0:*)~%")
      (format fout "      type(cint_env_vars), intent(in) :: envs~%")
      (format fout "      integer,  intent(in)    :: gout_empty~%")
      (format fout "      integer  :: nf, ix, iy, iz, n~%")
      (emit-g-decls fout ng)
      (when *f90-gptr* (emit-gptr-decls fout *f90-gptr*))
      (format fout "      real(dp) :: s(0:~a)~%" (1- (expt 3 tot-bits)))
      (dri-decl-f90 fout bra-i "i")
      (dri-decl-f90 fout (append op ket-j) "j")
      (dri-decl-f90 fout bra-k "k")
      (giao-decl-f90 fout bra-i (append op ket-j))
      (format fout "~%")
      (dri-init-f90 fout bra-i "i")
      (dri-init-f90 fout (append op ket-j) "j")
      (dri-init-f90 fout bra-k "k")
      (giao-init-f90 fout bra-i (append op ket-j))
      (format fout "      nf = envs%nf~%")
      (format fout "      g0 = 0~%")
      (loop for i in (range ng) do
            (format fout "      g~a = g~a + envs%g_size * 3~%" (1+ i) i))

      (dump-combo-braket fout (f90-fmt-1e-k-fn
                                (format nil "envs%i_l+~a" i-len)
                                (format nil "envs%j_l+~a" (+ op-len j-len)))
                         (lambda (&rest a) (declare (ignore a))
                           (error "no operator on the k side"))
                         (lambda (&rest a) (declare (ignore a))
                           (error "no l index in a three-centre 1e integral"))
                         k-rev '() '() 0)
      (dump-combo-braket fout (f90-fmt-1e-i3-fn) (f90-fmt-1e-op3-fn) (f90-fmt-1e-j3-fn)
                         i-rev op-rev j-rev k-len)

      (when *f90-gptr* (emit-gptr-binds fout *f90-gptr*))
      (format fout "~%      do n = 0, nf - 1~%")
      (format fout "         ix = idx(0+n*3)~%")
      (format fout "         iy = idx(1+n*3)~%")
      (format fout "         iz = idx(2+n*3)~%")
      (dump-s-for-nroots-f90 fout tot-bits 1)
      (gen-f90-block-with-empty fout flat-script)
      (format fout "      end do~%")
      (format fout "   end subroutine CINTgout1e_~a~%~%" intname)
      (setf *f90-gptr* nil))))

(defun gen-f90-int3c1e (fout intname raw-infix)
  (destructuring-bind (op bra-i ket-j bra-k ket-l)
    (split-int-expression raw-infix)
    (declare (ignore ket-l))
    (let* ((i-rev (effect-keys bra-i))
           (j-rev (reverse (effect-keys ket-j)))
           (k-rev (effect-keys bra-k))
           (op-rev (reverse (effect-keys op)))
           (op-len (length op-rev))
           (i-len (length i-rev))
           (j-len (length j-rev))
           (k-len (length k-rev))
           (tot-bits (+ i-len j-len op-len k-len))
           (raw-script (eval-int raw-infix))
           (flat-script (flatten-raw-script (last1 raw-script)))
           (sf (cadr raw-script))
           (goutinc (length flat-script))
           (e1comps (if (eql sf 'sf) 1 4))
           (tensors (if (eql sf 'sf) goutinc (/ goutinc 4)))
           (int1e-type (cond ((member 'nuc raw-infix) "INT3C1E_NUC")
                             ((or (member 'rinv raw-infix)
                                  (member 'nabla-rinv raw-infix)) "INT3C1E_RINV")
                             (t "INT3C1E_OVLP")))
           (ngdef (format nil "[~d, ~d, ~d, 0, ~d, ~d, 0, ~d]"
                          i-len (+ op-len j-len) k-len tot-bits e1comps tensors))
           (fac (factor-of raw-infix)))

      (format fout "   ! <(~{~a ~}i) (~{~a ~}~{~a ~}j) (~{~a ~}k)>~%" bra-i op ket-j bra-k)
      (gen-f90-optimizer fout intname ngdef "cint_all_3c1e_optimizer")
      (gen-f90-gout3c1e fout intname raw-infix flat-script)
      (dolist (form '(("cart" . "C2S_CART_3C2E1") ("sph" . "C2S_SPH_3C2E1")))
        (format fout "   function ~a_~a(out, dims, shls, atm, natm, bas, nbas, env, ws) &~%"
                intname (car form))
        (format fout "         result(has_value)~%")
        (format fout "      real(dp), intent(inout) :: out(0:)~%")
        (format fout "      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas~%")
        (format fout "      integer,  target        :: atm(0:), bas(0:)~%")
        (format fout "      real(dp), target        :: env(0:)~%")
        (format fout "      type(cint_ws), intent(inout) :: ws~%")
        (format fout "      logical :: has_value~%")
        (format fout "      type(cint_env_vars) :: envs~%")
        (format fout "      integer, parameter :: ng(0:7) = ~a~%~%" ngdef)
        (format fout "      call cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)~%")
        (format fout "      envs%f_gout => CINTgout1e_~a~%" intname)
        (unless (= fac 1)
          (format fout "      envs%common_factor = envs%common_factor * (~a_dp)~%" fac))
        (format fout "      has_value = cint_3c1e_drv(out, dims, envs, ws, ~a, ~a)~%"
                (cdr form) int1e-type)
        (format fout "   end function ~a_~a~%~%" intname (car form))))))

;;; ============================================================== the module
;;;
;;; One entry point for every arity, dispatching exactly as gen-cint does and
;;; on the same predicates.  auto_intor.cl is fed to this unchanged except for
;;; the call name and the output file, which is the point: one description
;;; list, two languages.

(defun f90-arity (raw)
  (cond ((int3c1e? raw)      '3c1e)
        ((int1e-grids? raw)  'grids)
        ((one-electron-int? raw) '1e)
        ((int4c2e? raw)      '4c2e)
        ((int3c2e? raw)      '3c2e)
        ((int2c2e? raw)      '2c2e)
        (t (error "unclassifiable integral"))))

(defparameter *f90-skipped* '())

(defun gen-f90-cint (filename modname &rest items)
  "Write one Fortran module for one of auto_intor.cl's output files."
  (let ((kept (remove-if (lambda (it)
                           (let ((a (f90-arity (cadr it))))
                             (declare (ignore a))
                             nil))
                         items)))
    (with-open-file (fout (mkstr filename) :direction :output :if-exists :supersede)
      (format fout "!~%! Generated by scripts/f90-backend.cl from the integral~%")
      (format fout "! descriptions in scripts/auto_intor.cl -- do not edit.~%!~%")
      (format fout "! The symbolic layer that produced these is the same one that~%")
      (format fout "! produces libcint's C: only the emitter differs.  The _spinor~%")
      (format fout "! entry points come out of the same description as _cart and~%")
      (format fout "! _sph, which is the whole argument for owning the generator.~%!~%")
      (format fout "module ~a~%" modname)
      (format fout "   use cint_const,     only: dp~%")
      (format fout "   use cint_envs~%")
      (format fout "   use cint_workspace, only: cint_ws~%")
      (format fout "   use cint_g1e~%")
      (format fout "   use cint_g2e~%")
      (format fout "   use cint_g2e_ops~%")
      (format fout "   use cint_1e,        only: cint_1e_drv~%")
      (format fout "   use cint_1e_spinor, only: cint_1e_spinor_drv~%")
      (format fout "   use cint_2e,        only: cint_2e_drv~%")
      (format fout "   use cint_2e_spinor, only: cint_2e_spinor_drv~%")
      (format fout "   use cint_3c2e,      only: cint_3c2e_drv, cint_2c2e_drv~%")
      (format fout "   use cint_3c2e_spinor, only: cint_3c2e_spinor_drv~%")
    ;; Each arity's optimizer builder lives in the module that owns its envs
    ;; init, because cint_opt has to stay underneath all of them.
    (format fout "   use cint_opt,       only: cint_del_optimizer~%")
    ;; cint_g1e and cint_g2e are already used bare above, so naming their
    ;; optimizer builders again in an ONLY clause adds no visibility -- the
    ;; bare use has already brought in everything those modules export.  ifx
    ;; says so out loud, remark #6536, once per module per generated file:
    ;; seventeen files, two modules each, thirty-four lines of noise in every
    ;; Intel build.  The builders stay reachable through the bare use; only
    ;; the redundant second mention goes.
    (format fout "   use cint_3c1e,      only: cint_all_3c1e_optimizer~%")
    (format fout "   use cint_1e_grids,  only: cint_all_1e_grids_optimizer~%")
    (format fout "   use cint_1e_grids_spinor, only: cint_1e_grids_spinor_drv~%")
    (format fout "   use cint_1e_grids,  only: cint_1e_grids_drv, GRID_BLKSIZE, &~%")
    (format fout "                             cint_init_int1e_grids_envvars, &~%")
    (format fout "                             cint_nabla1i_grids, cint_nabla1j_grids, &~%")
    (format fout "                             cint_x1i_grids, cint_x1j_grids~%")
    (format fout "   use cint_3c1e,      only: cint_3c1e_drv, cint_init_int3c1e_envvars, &~%")
    (format fout "                             INT3C1E_OVLP, INT3C1E_RINV, INT3C1E_NUC~%")
      (format fout "   implicit none~%")
      (format fout "   private~%~%")
      (dolist (item kept)
        (let ((nm (mkstr (car item))))
          ;; No _spinor for the derivative 2c2e forms or for any 3c1e form.
          ;; Both refusals are upstream's: CINT3c1e_spinor_drv is a stub that
          ;; exits, and the C's own generator emits a 2c2e spinor entry point
          ;; that prints "&c2s_sf_1e_spinor not implemented" and returns zero.
          ;; The plain int2c2e_spinor is real and lives in cint_3c2e_spinor.
          (if (member (f90-arity (cadr item)) '(2c2e 3c1e))
              (format fout "   public :: ~a_cart, ~a_sph~%" nm nm)
              (format fout "   public :: ~a_cart, ~a_sph, ~a_spinor~%" nm nm nm))
          ;; every integral gets an optimizer, as in the C, whatever its arity
          (format fout "   public :: ~a_optimizer~%" nm)
          ;; The gout, but only for the four the F12 entry points reuse by
          ;; name.  Exporting all 197 of them costs 2% on the two-electron
          ;; path -- measured -- because it denies LTO the chance to treat
          ;; them as local.  The C gives them all external linkage and pays
          ;; nothing, having no LTO to lose.
          (when (member nm '("int2e_ip1" "int2e_ipip1" "int2e_ipvip1" "int2e_ip1ip2")
                        :test #'string=)
            (format fout "   public :: CINTgout2e_~a~%" nm))))
      (format fout "~%contains~%~%")
      (dolist (item kept)
        (let ((raw (cadr item)) (nm (mkstr (car item))))
          (case (f90-arity raw)
            (3c1e (gen-f90-int3c1e fout nm raw))
            (grids (gen-f90-int1e-grids fout nm raw))
            (1e   (gen-f90-int1e   fout nm raw))
            (4c2e (gen-f90-int4c2e fout nm raw))
            (3c2e (gen-f90-int3c2e fout nm raw))
            (2c2e (gen-f90-int2c2e fout nm raw)))))
      (format fout "end module ~a~%" modname))))

(defun f90-report-skipped ()
  (format t "~%skipped (runtime not ported yet):~%")
  (dolist (e (reverse *f90-skipped*))
    (format t "  ~a  ~a  (~a)~%" (first e) (second e) (third e)))
  (format t "  ~a total~%" (length *f90-skipped*)))
