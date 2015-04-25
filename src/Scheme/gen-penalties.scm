(load "matrix.scm")
(load "~/lib/sets.scm")

; quicky to generate penalty table for lpa-align
;
(define amb_defs '( (A  (A))      (C  (C))    (G  (G))       (T  (T))
                  (M  (A  C))   (R  (A  G)) (W  (A  T))    (S  (C  G))
                  (Y  (C  T))   (K  (G  T)) (V  (A  C  G)) (H  (A  C  T))
                  (D  (A  G  T)) (B  (C  G  T)) (N  (A  C  G  T))
                 ))

(define codes (map car amb_defs))
(define num_codes (length codes))

(define is-amb?
  (lambda (x)
     (not (or (equal? x 'A) (equal? x 'C) (equal? x 'G) (equal? x 'T)))))

(define ambig_penalty 0)
(define subst_penalty 4)
(define undef_penalty -1)  ; should never show up

(define sym->char
  (lambda (sym)
    (char-upcase (string-ref (symbol->string sym) 0))))

(define match-code
  (lambda (x y)
    (if (null? (inter (cadr (assoc x amb_defs)) (cadr (assoc y amb_defs))))
        subst_penalty
        (if (or (is-amb? x) (is-amb? y))
            ambig_penalty
            0))))

(define gen-penalty-matrix
  (lambda ()
    (let ( (mat (make-matrix 256 256 undef_penalty)) )
      (for-each 
        (lambda (i)
            (for-each 
              (lambda (j)
                (let ( (val (match-code i j))
                       (a (sym->char i))  ; uppercase
                       (b (sym->char j))  ; also uppercase
                       )
                  (let (
                        (ua (char->integer a))
                        (ub (char->integer b))
                        (da (char->integer (char-downcase a)))
                        (db (char->integer (char-downcase b)))
                        )
                    (set-matrix! mat ua ub val)
                    (set-matrix! mat ua db val)
                    (set-matrix! mat da ub val)
                    (set-matrix! mat da db val)
                    )))
              codes)
              )
        codes)
      mat
      )))

(define penalty_matrix (gen-penalty-matrix))

;
; (penalty-value a b) - returns match score for characters a & b
;

(define penalty-value
  (lambda (a b)
    (matrix-ref penalty_matrix (char->integer a) (char->integer b))))


(define lets (map sym->char codes))



