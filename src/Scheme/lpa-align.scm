(load "~/lib/sets.scm")
(load "prqueue.scm")
(load "matrix.scm")

;
; Constants
;
(define indel_penalty 5)    ; penalty for insertion or deletion
(define subst_penalty 4)    ; penalty for substitution
(define ambig_penalty 0)    ; penalty for ambiguity substition
(define undef_penalty -1)   ; marks undef locs in penalty matrix

;(define up_pnt  1)          ; directions for backpnt matrix
;(define diag_pnt 2)
;(define left_pnt 3)

(define indel_ind #\-)      ; printable values used in mid array
(define ambig_ind #\|)
(define wrong_ind #\*)
(define perft_ind #\space)
(define indel_spc #\space)  ; for indels in top and bot arrays

(define unseen 100000000)   ; marks unseen locations in dynamic matrix d.
                            ; this should be maximally huge and positive

(define penalty_matrix '())

;(lpa-align "seq1" "seq2") returns an alignment matrix (dsynamic programming)
;                          seq1 & seq2 can be upper or lower case; casefolding
;                          is taken care of in the penalty matrix
;
; Note to self: this uses a global for the priority queue - maybe it should
;               be a local, to allow multiple calls 
;
(define lpa-align
  (lambda (seq1 seq2)
    (if (null? penalty_matrix)
        (set! penalty_matrix (gen-penalty-matrix)))    
    (let ( (m (string-length seq1)) (n (string-length seq2)) )
      (let ( (d                (make-matrix (+ m 1) (+ n 1) (- unseen)))
             (backpnt          (make-matrix (+ m 1) (+ n 1) #f))
             (max_right_offset (estim-max-right-offset n m))
             (max_down_offset  (estim-max-down-offset n m))
             )
        (pq-init)
        (pq-update unseen '(0 0))
        (let fill ( (ij (pq-remove)) )
          (if (not ij) (error "lpa-align: priority queue prematurely empty"))
          (let ( (i (caadr ij)) (j (cadadr ij)) )
          (disp "pq-remove " ij  "i = " i "j = " j)
            (if (or (< i m) (< j n))
                (begin
                  (set-matrix! d i j (- (matrix-ref d i j)))
                  (if (= (matrix-ref d i j) unseen) 
                      (set-matrix! d i j 0))
                  ;
                  ; right - [i,j+1]
                  ;
                  (if (and (< j n) (< (matrix-ref d i (+ j 1)) 0))
                      (let ( (dist (+ (matrix-ref d i j)
                                      (cond ((= i m)                      0)
                                            ((and (= i 0) 
                                                  (< j max_right_offset)) 0)
                                            (else            indel_penalty)))))
                        (if (pq-update dist (list i (+ j 1)))
                            (begin
                              (set-matrix! d i (+ j 1) (- dist))
                              (set-matrix! backpnt i (+ j 1) 'left_pnt)
                              (disp "right" dist (matrix-ref d i (+ j 1)))
                              ))
                        ))
                  ;
                  ; down - [i+1,j]
                  ;
                  (if (and (< i m) (< (matrix-ref d (+ i 1) j) 0))
                      (let ( (dist (+ (matrix-ref d i j)
                                      (cond ((= j n)                      0)
                                            ((and (= j 0) 
                                                  (< i max_down_offset))  0)
                                            (else            indel_penalty)))))
                        (if (pq-update dist (list (+ i 1) j))
                            (begin
                              (set-matrix! d (+ i 1) j (- dist))
                              (set-matrix! backpnt (+ i 1) j 'up_pnt)
                              (disp "down" dist (matrix-ref d (+ i 1) j))
                                      
                              ))
                        ))
                  ;
                  ; diag - [i+1,j+1]
                  ;
                  (if (and (< i m) (< j n) (< (matrix-ref d (+ i 1) (+ j 1)) 0))
                      (let ( (dist (+ (matrix-ref d i j)
                                      (penalty-value (string-ref seq1 i)
                                                     (string-ref seq2 j)))) )
                        (if (pq-update dist (list (+ i 1) (+ j 1)))
                            (begin
                              (set-matrix! d (+ i 1) (+ j 1) (- dist))
                              (set-matrix! backpnt (+ i 1) (+ j 1) 'diag_pnt)
                              (disp "diag" dist (matrix-ref d (+ i 1) (+ j 1)))
                              ))
                        ))
                  
                  ; recurse - get next
                  (fill (pq-remove))
                  )
                (begin
                  (dump-matrix d)
                  (disp "backpnt")
                  (dump-matrix backpnt)
                  (gen-alignment backpnt seq1 seq2)
                  )
                ))
          )))))

;
; gen-alignment takes the 2-d backpnt matrix and works backwards to
;               produce the top, mid and bot arrays
;
(define gen-alignment
  (lambda (backpnt seq1 seq2)
    (disp "generating alignment")
    (let backp ( (i       (string-length seq1)) 
                 (j       (string-length seq2)) 
                 (algnmt '() ) )
      (disp 'backp i j algnmt)
      (if (or (> i 0) (> j 0))
          (cond ( (eq? (matrix-ref backpnt i j) 'up_pnt)

                  (backp (- i 1) j
                         (cons (list (string-ref seq1 (- i 1))
                                     indel_ind indel_spc)
                               algnmt)) )
                
                ( (eq? (matrix-ref backpnt i j) 'left_pnt)

                  (backp i (- j 1)
                         (cons (list indel_spc indel_ind
                                     (string-ref seq2 (- j 1)))
                               algnmt)) )

                ( (eq? (matrix-ref backpnt i j) 'diag_pnt)

                  (let ( (x (string-ref seq1 (- i 1)))
                         (y (string-ref seq2 (- j 1))) )
                    (backp (- i 1) (- j 1)
                         (cons (list x
                                     (cond ( (char=? x y)  perft_ind )
                                           ( (is-amb? x)   ambig_ind )
                                           ( (is-amb? y)   ambig_ind )
                                           ( else          wrong_ind ) )
                                     y)
                               algnmt))) )
                (else (disp "backp error: " i j (matrix-ref backpnt i j)))
            )
          (begin
            (disp "done generating alignment")
            algnmt)
          ))))

; these things are a kludge to prevent a complete disconnect, until I can
; come up with something better.
;
(define estim-max-right-offset
  (lambda (n m)
    (if (> n m) (- n m) (floor (/ n 10)) )))

(define estim-max-down-offset
  (lambda (n m)
    (if (> m n) (- m n) (floor (/ m 10)) )))


;
;  Penalty matrix 
;
;(conv-chars s) - returns structure s with all symbols -> chars
;
(define conv-chars
  (lambda (s)
    (cond ( (null? s)    '() )
          ( (symbol? s)  (char-upcase (string-ref (symbol->string s) 0)) )
          ( (pair? s)    (cons (conv-chars (car s)) (conv-chars (cdr s))) )
          ( else         (error "bad thing occured in conv-chars") )
          )))


;
; quicky to generate penalty table for lpa-align
;
(define amb_defs
            (conv-chars
               '( (A  (A))      (C  (C))    (G  (G))       (T  (T))
                  (M  (A  C))   (R  (A  G)) (W  (A  T))    (S  (C  G))
                  (Y  (C  T))   (K  (G  T)) (V  (A  C  G)) (H  (A  C  T))
                  (D  (A  G  T)) (B  (C  G  T)) (N  (A  C  G  T))
                 )))

(define codes (map car amb_defs))

(define is-amb?
  (lambda (x)
    (let ( (c (char-upcase x)) )
      (not (or (equal? c #\A) (equal? c #\C) (equal? c #\G) (equal? c #\T))))))

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
        (lambda (a)
            (for-each 
              (lambda (b)
                (let ( (val (match-code a b))
                       (ua (char->integer a))  ; a and b are uppercase
                       (ub (char->integer b))
                       (da (char->integer (char-downcase a)))
                       (db (char->integer (char-downcase b)))
                        )
                    (set-matrix! mat ua ub val)
                    (set-matrix! mat ua db val)
                    (set-matrix! mat da ub val)
                    (set-matrix! mat da db val)
                    ))
              codes)
              )
        codes)
      mat
      )))

;
; (penalty-value a b) - returns match score for characters a & b
;

(define penalty-value
  (lambda (a b)
    (matrix-ref penalty_matrix (char->integer a) (char->integer b))))

;
; (string-upcase str) returns a new string with all lowercase chars changed
;                     to uppercase
(define string-upcase 
  (lambda (str)
    (let ( (n (string-length str)) )
      (let ( (upstr (make-string n)) )
        (let cp ( (i 0) )
          (if (>= i n)
              upstr
              (begin 
                (string-set! upstr i (char-upcase (string-ref str i)))
                (cp (+ i 1))
                )))))))

(define disp
  (lambda lst
    (let dspl ((l lst))
      (if (null? l)
          (newline)
          (begin
            (display (car l))
            (display " ")
            (dspl (cdr l))
            )))))

