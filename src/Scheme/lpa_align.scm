;
; Priority queue, implemented as a binary tree for now
;
;(pq-update key val)  - inserts key and val
;(pq-remove)  -> returns (key val) or #f

(define pq_root '())

(define pq-update
  (lambda (key val)
    (if (null? pq_root)
        (begin 
          (set! pq_root (make-pq-node key '() '() val))
          #t)
        (pq-insert pq_root key val))))

(define pq-remove
  (lambda ()
    (if (null? pq_root)
        #f
        (let ( (r (pq-lremove pq_root)) )
          (set! pq_root (car r))
          (cdr r)))))



; (pq-insert node key val)
;
(define pq-insert
  (lambda (node key val)
    ;(write (list 'pq-insert node key val)) (newline)
    (if (< key (pq-key node))
        (if (null? (pq-left node))
            (begin 
               (set-pq-left! node (make-pq-node key '() '() val))
               #t)
            (pq-insert (pq-left node) key val))
        (if (null? (pq-right node))
            (begin
              (set-pq-right! node (make-pq-node key '() '() val))
              #f)
            (begin
              (pq-insert (pq-right node) key val)
              #f))
            
        )))

;
; (pq-lremove node) returns (node key val))  where (key val) are from the
;                   leftmost node removed, and node points to the result
;                   binary tree
;
(define pq-lremove
  (lambda (node)
    (if (null? (pq-left node))
        (list (pq-right node) (pq-key node) (pq-val node))
        (let ( (r (pq-lremove (pq-left node))) )
          ;(write (list 'r 'is r)) (newline)
          (set-pq-left! node (car r))
          (cons node (cdr r))
          ))))
;
;(pq-dump node)
;
(define pq-dump
  (lambda (node)
    (if (not (null? node))
        (begin
          (pq-dump (pq-left node))
          (write (list (pq-key node) (pq-val node))) (newline)
          (pq-dump (pq-right node))
          ))))


;
;(make-pq-node key left right val)
;
(define make-pq-node 
  (lambda (key left right val)
    (list key (cons left right) val)))

(define pq-key (lambda (node) (car node)))

(define pq-left (lambda (node) (caadr node)))

(define pq-right (lambda (node) (cdadr node)))

(define pq-val (lambda (node) (caddr node)))

(define set-pq-left!
  (lambda (node left)
    (set-car! (cadr node) left)))

(define set-pq-right!
  (lambda (node left)
    (set-cdr! (cadr node) left)))

; pq checks
;(define x (make-pq-node 123 '() '() 'high))
(pq-update 69 'p1)
(pq-update 73 'p2)
(pq-update 54 'p3)
(pq-update 141 'p5)
(pq-update 100 'p4)

(pq-dump pq_root)

