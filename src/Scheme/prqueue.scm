;
; Priority queue, implemented as a binary tree for for lpa-align.  This
; is kind of weird because it also maintains a hash to allow lookups by 
; val as well.
;
; Dependencies:  mzscheme hash tables
;

;
;(pq-init)
;(pq-update key i j)  - inserts key and (i j) - returns #t if this key is
;                       is the new minimum, otherwise it returns false
;(pq-remove)  -> returns (key val) or #f

(define pq_root     '())          ; binary tree root
(define pq_val_hash '())          ; hash table with [i,j] value as key
(define pq_hash_off  1000000)

;
; (pq-init) initializes (resets) the priorty queue data structures.
;
(define pq-init 
  (lambda () 
    (set! pq_root '())
    (set! pq_val_hash (make-hash-table))
    ))

(define pq-update
  (lambda (key i j)
    (let ( (hkey (pq-hash-key i j)) )
      (if (null? pq_root)
          (begin
            (set! pq_root (make-pq-node key '() '() (i j)))
            (hash-table-put! pq_val_hash hkey pq_root)
            )
          (let ( (p (hash-table-get pq_val_hash hkey (lambda () #f)) ) )
            (if p
                ; found it - if key in pq tree is < key, do nothing, return #f)
                ; otherwise, update key in tree, peform sort
                (begin
                  (display "exists: ")
                  (if (< (pq-key p) key)
                      (begin (display "not updated") (newline)
                             #f)
                      (begin
                        (set-pq-key! p key)
                        ; resort
                        #t))
                  )
                
      

;
; (pq-hash-key i j) returns a single hash key value out of i j
;
(define pq-hash-key
  (lambda (i j)
    (+ (* pq_hash_off i) j)))

  
(define old-pq-update
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

(define set-pq-key!
  (lambda (node key)
    (set-car! node key)))

; pq checks
;(define x (make-pq-node 123 '() '() 'high))
;(pq-update 69 'p1)
;(pq-update 73 'p2)
;(pq-update 54 'p3)
;(pq-update 141 'p5)
;(pq-update 100 'p4)

;(pq-dump pq_root)

