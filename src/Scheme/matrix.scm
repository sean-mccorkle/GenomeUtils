;  matrix.scm - matrix routines
;
; (make-matrix n m init_val)    returns a matrix
; (matrix-ref mat i j)          returns value at i,j
; (set-matrix! mat i j val)     sets value at i,j
; (dump-matrix mat)             prints out matrix

(define make-matrix
  (lambda (n m init_val)
    (list n m (make-vector (* n m) init_val))))

(define matrix-ref 
  (lambda (mat i j)
    (if (and (>= i 0) (< i (num-rows mat)) (>= j 0) (< j (num-cols mat)))
        (vector-ref (caddr mat) (+ (* i (num-cols mat)) j))
        (error "matrix ref out of range")
        )))

(define set-matrix! 
  (lambda (mat i j val)
    (if (and (>= i 0) (< i (num-rows mat)) (>= j 0) (< j (num-cols mat)))
        (vector-set! (caddr mat) (+ (* i (num-cols mat)) j) val)
        (error "matrix ref out of range")
        )))
   
; (num-rows matrix) - returns number of rows
; (num-cols matrix) - returns number of columns

(define num-rows  (lambda (mat) (car mat)) )
(define num-cols  (lambda (mat) (cadr mat)) )

(define dump-matrix
  (lambda (mat)
    (let ( (n (num-rows mat)) (m (num-cols mat)) )
      (let l1 ( (i 0) )
        (if (< i n)
            (begin
              (let l2 ( (j 0) )
                (if (< j m)
                    (begin
                      (display (matrix-ref mat i j))
                      (display " ")
                      (l2 (+ j 1))
                      )
                    (newline)))
              (l1 (+ i 1))
              ))))))
     
;(define mx (make-matrix 4 5 2))
;(dump-matrix mx)
;(newline)

;(set-matrix! mx 0 0 7)
;(dump-matrix mx)
;(newline)

;(set-matrix! mx 0 4 0)
;(set-matrix! mx 0 4 1)
;(dump-matrix mx)
;(newline)
