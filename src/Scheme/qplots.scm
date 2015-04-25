; Program:      qplots.scm
; Programmer:   Sean R. McCorkle
; Language:     Scheme (MZ Scheme)
;
; Description:  Reads phred quality scores (fasta-form) and generates 
;               small miniature plots of q vs position or histograms of
;               q for each sequence.  Q values less than the standard
;               threshold (20, corresponding to 1% error probability)
;               are printed in grey, values exceeding the threshold are
;               printed in black.
;               Currently, this generates postscript output on stdout.
;
; Usage:        mzscheme -r qplots.scm [-h] [-H]   qualfile
;
; $Id: qplots.scm,v 0.2 2007/01/21 18:38:46 mccorkle Exp mccorkle $
;
(require (lib "cmdline.ss"))

(load "~/lib/io.scm")
(load "~/lib/plot_ev2.scm")
(load "~/lib/hist.scm")


(define quals '())
(define quals_filename "")
(define plot_hists  #f)

(define num_rows 10)                ; 10 rows of plots per page
(define num_cols 4)                 ; four columns of plots across a page
(define q_max 70) ; 60              ; maximum value of q expected
(define q_thresh 20)                ; 
(define q_frame_width 120)
(define q_frame_height 40)
;(define q_frame_width 500)
;(define q_frame_height 200)

(define left_margin 20)
(define top_margin 680)
(define horz_space 20)
(define vert_space 25)

;(define top_margin 500)

(define max_len 1000)  ; 200

(define truncate-at-blank
  (lambda (str)
    (let ( (n (string-length str)) )
      (let find ( (i 0) )
        (if (and (< i n) (not (char=? (string-ref str i) #\space)))
            (find (+ i 1))
            (substring str 0 i)
            )))))

(define ro-vector
  (lambda (n)
    (let ( (v (make-vector n)) )
      (let fill ( (i 0) )
        (if (< i n)
            (begin
              (vector-set! v i i)
              (fill (+ i 1))
              )
            v)))))

(define truncate-vector
  (lambda (v n)
    (if (<= (vector-length v) n)
        v
        (let ( (new (make-vector n)) )
          (let fill ( (i 0) )
            (if (< i n)
                (begin
                  (vector-set! new i (vector-ref v i))
                  (fill (+ i 1))
                  )
                new))))))
;
; (next-pos vec start n pred? returns the first index i of vec, after start, 
;   where (pred? vec[i]) is #t
;   
(define next-pos
  (lambda (vec start n pred?)
    (let srch ( (i start) )
      (if (< i n)
          (if (pred? (vector-ref vec i))
              i
              (srch (+ i 1)))
          i))))

; returns (start stop) or #f
;
(define next-interval
  (lambda (vec start n pred?)
    (let ( (a (next-pos vec start n pred?)) )
      (if (< a n)
          (list a (next-pos vec a n (lambda (v) (not (pred? v)))))
          #f)
      )))

(define subvector
  (lambda (v a b)
    (let ( (n (vector-length v)) )
      (if (< a n)
          (let ( (c (if (< b n) b n)) )
            (let ( (subv (make-vector (- c a))) )
              (let fill ( (i 0) (j a) )
                (if (< j c)
                    (begin
                      (vector-set! subv i (vector-ref v j))
                      (fill (+ i 1) (+ j 1)))
                    subv))))
          (make-vector 0)
          ))))

  
;(define v #( 1 0 1 0 1 0 5 6 5 6 6 5))

;(next-pos v 0 12 (lambda (x) (>= x 5)))

;(next-pos v 0 12 (lambda (x) (> x 5)))
;(next-pos v 8 12 (lambda (x) (> x 5)))
;(next-pos v 10 12 (lambda (x) (<= x 5)))
;(next-pos v 5 12 (lambda (x) (= x 100)))

;(next-interval v 0 12 (lambda (x) (> x 4)))
;(next-interval v 0 12 (lambda (x) (> x 5)))
;(next-interval v 8 12 (lambda (x) (> x 5)))
;(next-interval v 8 12 (lambda (x) (< x 5)))
;(next-interval v 0 12 (lambda (x) (< x 5)))

(define origin
  (lambda (i)
    (list (+ (* (+ q_frame_width horz_space) (modulo i num_cols)) left_margin)
          (- top_margin (* (+ q_frame_height vert_space) (floor (/ i num_cols))))
          )))
          

(define plot-intervals
  (lambda (fr qv xs n pred? graylevel)
    (set-gray graylevel)
    (let pl ( (int (next-interval qv 0 n pred?)) )
      (if int
          (let ( (a (car int)) (b (cadr int)) )
            (bar-graph fr a b (subvector qv a b))
            (pl (next-interval qv b n pred?))
            )
          ))))
           
(define plot-qs
  (lambda (qv n x_orig y_orig name)
    (let ( (xs (ro-vector n))
           (fr (create-frame x_orig (+ x_orig q_frame_width)
                             y_orig (+ y_orig q_frame_height)
                             (range 0 max_len)
                             (range 0 q_max))) 
             )
      ;(bar-graph fr 0 n qv)
      (plot-intervals fr qv xs n (lambda (x) (>= x 20)) 100)
      (plot-intervals fr qv xs n (lambda (x) (< x 20)) 50)
      (set-gray 30.0)
      (frame-line fr 'solid 0 q_thresh max_len q_thresh)
      (set-gray 100)
      (frame-line fr 'solid 0 0 max_len 0)
      (frame-line fr 'solid 0 0 0 q_max)

      (for-each 
            (lambda (p) (frame-rline fr 'solid p 0 0 -3))
            '(100 200 300 400 500 600 700 800 900) )

      (scalefont 8)
      (text (+ x_orig q_frame_width) (+ y_orig q_frame_height 8)
            'right 'top name)
      )))

(define plot-hist
  (lambda (qv n x_orig y_orig name)
    (let ( (xs (ro-vector q_max))
           (h  (hist1d (vector->list qv) 0 q_max q_max)) )
      (let ( (yrange (find-range h #f)) )
        (let ( (fr (create-frame x_orig (+ x_orig q_frame_width)
                             y_orig (+ y_orig q_frame_height)
                             (range 0 q_max) yrange)) )
          (set-gray 50.0)
          (bar-graph fr 0 q_thresh (subvector h 0 q_thresh))
          (set-gray 100.0)
          (bar-graph fr q_thresh q_max (subvector h q_thresh q_max))
;      (plot-intervals fr qv xs n (lambda (x) (>= x 20)) 100)
;      (plot-intervals fr qv xs n (lambda (x) (< x 20)) 50)
          (set-gray 30.0)
          ;(frame-line fr 'solid q_thresh 0 q_thresh (range-high yrange))
          (set-gray 100)
          (frame-line fr 'solid 0 0 q_max 0)
          (frame-line fr 'solid 0 0 0 (range-high yrange))
          (frame-rline fr 'solid 20 0 0 -3)
          (frame-rline fr 'solid 40 0 0 -3)
          (frame-rline fr 'solid 60 0 0 -3)
          (scalefont 8)
          (text (+ x_orig q_frame_width) (+ y_orig q_frame_height 8)
               'right 'top name)
          )))))


(define finish-page
  (lambda ()
    (disp "showpage")
    ))

(define new-page
  (lambda () (disp "% new page") ))
  
(define process
  (lambda (qs name counter)
    (if (>= counter 0)   ; don't bother if we haven't got data yet from main loop
      (let ( (qvec (truncate-vector (list->vector (reverse qs)) max_len))
             (page_counter (modulo counter (* num_rows num_cols)))
             )

        (let ( (n (vector-length qvec))
               (orig (origin page_counter))
               )
          (disp "%process length " n)
          (let ( (x_orig (car orig))
                 (y_orig (cadr orig)) )
            (if plot_hists
                (plot-hist qvec n x_orig y_orig name)
                (plot-qs qvec n x_orig y_orig name)
                )
            )
          (if (= page_counter (- (* num_rows num_cols) 1))
              (begin
                (finish-page)
                (new-page)
                ))
          )))))
                     
                                          ;;;;;;;;;;;;;;;;
                                          ; Main Program ;
                                          ;;;;;;;;;;;;;;;;

(command-line "qplots" (current-command-line-arguments)
              (once-each
                 (("-H" "--hist") "histograms" (set! plot_hists #t))
                 )
              (args (fname) (set! quals_filename fname))
              )

(define f (open-input-file quals_filename))
(define counter -1)
(define seq_name "")
(init-canvas 'ps 800 1200)

;
; Main read loop - scan input, for each sequence, put quality scores into quals,
; and then invoke process
;

(let readl ( (line (getline f)) )
  (cond ( (eof-object? line)                (process quals seq_name counter) )
        ( (char=? (string-ref line 0) #\>)  (process quals seq_name counter)
                                            (set! quals '())
                                            (set! seq_name 
                                                  (truncate-at-blank 
                                                    (substring line 1 
                                                             (string-length line))))
                                            (set! counter (+ counter 1))
                                            (readl (getline f)) )
        ( else
                                            (let ( (s (open-input-string line)) )
                                               (let rds ( (n (read s)) )
                                                 (if (not (eof-object? n))
                                                     (begin
                                                       (set! quals (cons n quals))
                                                       (rds (read s))
                                                       ))))
                                            (readl (getline f)) )
        ))

(close-canvas)
(close-input-port f)

