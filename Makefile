REMOTE_DIR = /var/www/html/feb-2024
MATLAB=/Applications/MATLAB_R2024a.app/bin/matlab

all: qmd slides_pdf matlab zipit deploy

qmd:
	for f in $$(find . -name \*.qmd); do \
		echo "Rendering $$f\n"; \
		quarto render $$f; \
	done

slides_pdf:
	for f in $$(find _output -name folien.html); do \
		g=$$(echo $$f | sed "s/\.html/.pdf/"); \
		decktape $$f $$g; \
	done

MLX_FILES=\
02-stabtragwerke/aufgaben/musterloesung/bofem/fachwerk.mlx \
02-stabtragwerke/aufgaben/musterloesung/bofem/rahmen.mlx \
03-stabilitaet/folien/matlab/fachwerkstuetze.mlx

matlab:
	for f in $(MLX_FILES); do \
		echo $$f; \
		wd=$$(pwd); \
		d=$$(dirname $$f); \
		o=$$(basename $$f); \
		h=$$(echo $$o | sed "s/\.mlx/.html/"); \
		od=$$wd/_output/$$d; \
		of=$$od/$$h; \
		cd $$d; \
		test -f $$of && rm $$of; \
		echo; \
		pwd; \
		$(MATLAB) -nodesktop -r "export('$$o', '$$of'); exit()"; \
		mkdir -p $$od; \
		cd $$wd; \
	done

zipit:
	for f in */folien/01-demos; do \
		cd $$f && bash zipit.sh; \
		mkdir -p ../../../_output/$$f; \
		mv *.zip ../../../_output/$$f; \
		cd -; \
	done

deploy:
	rsync -v -a _output maba@igc:$(REMOTE_DIR)

clean:
	rm -rf _output
