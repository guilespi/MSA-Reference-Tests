(ns multiple-seq-aligner.core
  (:require [pallet.stevedore :refer [script with-script-language]]
            [pallet.common.shell :as shell])
  (:require pallet.stevedore.bash))

(def basedir "/Users/guilespi/Documents/Development/interrupted/bioinformatics/")
(def balibase (str basedir "BAliBASE2"))
(def baliscore "")

(defmacro find-cmd
  [base-dir pattern]
  `(chain-and ("cd" ~base-dir) ("find" "-name" ~pattern)))

(defmacro exec-script
  [s]
  `(shell/bash 
    (with-script-language :pallet.stevedore.bash/bash
      (script
       ~s))))

(defn find-files
  [path pattern]
  (let [change-dir (str "cd " path)
        query (str "find" " -name " pattern)
        response (exec-script (chain-and (change-dir) (query)))]
    (if (= 0 (:exit response))
      (clojure.string/split (:out response) #"\n")
      (:err response))))

(defn fasta-map 
  []
  (println "Building fasta map")
  (let [path (str balibase "/fastas")
        file-list (find-files path "*")
        fastas (filter #(re-find #"\.india|\.tfa" %) file-list)]
    (reduce (fn [hm f] 
              (let [key (re-find #"^[^.]+" f)]
                (assoc hm key (str path "/" f)))) {} fastas)))

(defn reference-map 
  []
  (println "Building reference map")
  (let [path (str balibase "/references")
        file-list (find-files path "*")
        references (filter #(re-find #"\.msf" %) file-list)]
    (reduce (fn [hm f] 
              (let [key (re-find #"^[^.]+" f)]
                (assoc hm key (str path "/" f)))) {} references)))

(def dialign-path (str basedir "tools/dialign_package"))
(def output-dir "results")

(defn dialign
  [sequence-file file-name output-dir]
  (println "Building dialign alignments...")
  (let [program (str dialign-path "/src/dialign2-2")
        data-dir (str dialign-path "/dialign2_dir")
        aligner (format "%s -fn ./%s/%s_dialign -msf %s" program output-dir file-name sequence-file)
        command (format "export DIALIGN2_DIR=%s;%s" data-dir aligner)
        response (exec-script (command))]
    (when (= "" (:out response))
      (format "%s/%s_dialign.ms" output-dir file-name))))

(defn align-sequences
  [fastas algorithm]
  (reduce (fn [hm a] 
            (let [key (get a 0) 
                  file-name (clojure.string/replace key #"/" "_")
                  start-time (System/nanoTime)
                  msf-file (algorithm (get a 1) file-name output-dir)
                  stop-time (System/nanoTime)
                  response {:msf-file msf-file
                            :time (- stop-time start-time)}]
              (assoc hm key response))) {} fastas))

(def tcoffee-path (str basedir "tools/tcoffee/Version_9.03.r1318/bin/t_coffee"))

(defn clustal-to-msf
  [input-file output-file]
  (let [command (format "%s -other_pg seq_reformat -in %s -output msf" tcoffee-path input-file)
        response (exec-script (command))]
    (when (= 0 (:exit response))
      (let [content (:out response)]
        (spit output-file content)
        true))))

(defn tcoffee
  [sequence-file file-name output-dir]
  (println "Building tcoffee alignments...")
  (let [tool tcoffee-path
        output-file (format "./%s/%s_tcoffee" output-dir file-name)
        parameters (format " %s -output=msf -run_name=%s" sequence-file output-file)
        command (str tool parameters)
        response (exec-script (command))]
    (when (= 0 (:exit response))
      (str output-file ".msf"))))

(def poa-path (str basedir "tools/poaV2"))
(defn poa 
  [sequence-file file-name output-dir]
  (println "Building poa alignments...")
  (let [tool (str poa-path "/poa")
        output-file (format "./%s/%s_poa" output-dir file-name)
        matrix-file (str poa-path "/blosum80.mat")
        parameters (format " -read_fasta %s -clustal %s %s" sequence-file output-file matrix-file)
        command (str tool parameters)
        response (exec-script (command))]
    (when (= 0 (:exit response))
      (let [msf-file (str output-file ".msf")]
        (when (clustal-to-msf output-file msf-file)
          msf-file)))))

(def clustal-path (str basedir "tools/clustalw-2.1-macosx/clustalw2"))
(defn clustalw2 
  [sequence-file file-name output-dir]
  (println "Building clustalw2 alignments...")
  (let [tool clustal-path
        output-file (format "./%s/%s_clustalw2.msf" output-dir file-name)
        parameters (format " -INFILE=%s -ALIGN -QUIET -OUTFILE=%s -OUTPUT=GCG" sequence-file output-file)
        command (str tool parameters)
        response (exec-script (command))]
    (when (= 0 (:exit response))
      output-file)))

(def baliscorer-path (str basedir "tools/bali_score"))
(defn calculate-bali-score
  [reference alignment]
  (let [command (format "%s %s %s" 
                        baliscorer-path
                        reference
                        (:msf-file alignment))
        response (exec-script (command))]
    (when (= 0 (:exit response))
      (let [output (:out response)
            [_ sp-score] (re-find #"SP score= ([.0-9]+)" output)
            [_ tc-score] (re-find #"TC score= ([.0-9]+)" output)]
        {:sp-score sp-score
         :tc-score tc-score
         :time (/ (:time alignment) 1000000.0)}))))


(defn calculate-results
  [fastas references algs]
  (println "Exporting results to csv")
  (reduce (fn [hm f] 
            (let [key (get f 0)
                  refer (get references key)
                  row (reduce #(assoc %1 
                                 (get %2 0)
                                 (calculate-bali-score refer (get (get %2 1) key))) {} algs)]
              (assoc hm key row))) {} fastas))

(def export-keys [:referid 
                  :testid
                  :seqid 
                  :tcoffee-sp
                  :tcoffee-tc
                  :tcoffee-time
                  :poa-sp
                  :poa-tc
                  :poa-time
                  :clustalw2-sp
                  :clustalw2-tc
                  :clustalw2-time
                  :dialign-sp
                  :dialign-tc
                  :dialign-time])

(defn save-results
  [results]
  (with-open [wrtr (clojure.java.io/writer "./results/scores.csv")]
    (.write wrtr (str (clojure.string/join "," (map #(str "\"" (name %) "\"") export-keys)) "\n"))
    (doseq [r results]
      (let [key (get r 0)
            scores (get r 1)
            [referid testid seqid] (clojure.string/split key #"/")
            row {:tcoffee-sp (-> scores :tcoffee :sp-score)
                 :tcoffee-tc (-> scores :tcoffee :tc-score)
                 :tcoffee-time (-> scores :tcoffee :time)
                 :poa-sp (-> scores :poa :sp-score)
                 :poa-tc (-> scores :poa :tc-score)
                 :poa-time (-> scores :poa :time)
                 :clustalw2-sp (-> scores :clustalw2 :sp-score)
                 :clustalw2-tc (-> scores :clustalw2 :tc-score)
                 :clustalw2-time (-> scores :clustalw2 :time)
                 :dialign-sp (-> scores :dialign :sp-score)
                 :dialign-tc (-> scores :dialign :tc-score)
                 :dialign-time (-> scores :dialign :time)
                 :referid referid
                 :testid testid
                 :seqid seqid}
            values (reduce (fn[v k] (conj v (str "\"" (get row k) "\""))) [] export-keys)
            csv-row (str  (clojure.string/join "," values) "\n")]
          (.write wrtr csv-row)))))

(defn alignments
  []
  (let [fastas (fasta-map)
        dialign-res (align-sequences fastas dialign)
        clustalw-res (align-sequences fastas clustalw2)
        poa-res (align-sequences fastas poa)
        tcoffee-res (align-sequences fastas tcoffee)
        results (calculate-results fastas
                                   (reference-map)
                                   {:dialign dialign-res
                                    :clustalw2 clustalw-res
                                    :poa poa-res
                                    :tcoffee tcoffee-res})]
    (save-results results)))
