.PHONY: train perspective

SUBDIRS = train perspective

bin_SCRIPTS = \
	analyzeLog \
	computeFinalResult \
	cvplayer
CLEANFILES = $(bin_SCRIPTS)

analyzeLog: groundtruth/analyzeLog.py
	cp groundtruth/analyzeLog.py analyzeLog
	chmod +x analyzeLog

computeFinalResult: groundtruth/computeFinalResult.py
	cp groundtruth/computeFinalResult.py computeFinalResult
	chmod +x computeFinalResult

cvplayer: groundtruth/cvplayer.py
	cp groundtruth/cvplayer.py cvplayer
	chmod +x cvplayer

clean:
	rm -f   ../../training/icdar/.run_training.stamp \
		../../training/icdar/negatives.txt \
		../../training/icdar/positives.txt
	make -C train $@
	make -C perspective $@

