tests = tests
package = atropos
#pytestopts = -s -vv --show-capture=all
repo = jdidion/$(package)
desc = Release $(version)

all: clean install install_extra_requirements install_test_requirements test test_release_setup

clean:
	rm -Rf __pycache__
	rm -Rf .pytest_cache
	rm -Rf .eggs
	rm -Rf **/__pycache__/*
	rm -Rf **/*.c
	rm -Rf **/*.so
	rm -Rf **/*.pyc
	rm -Rf dist
	rm -Rf build
	rm -Rf atropos.egg-info
	rm -f .adapters
	rm -f .coverage
	rm -f coverage.xml
	rm -f MANIFEST

build:
	python setup.py build_ext -i
	python setup.py sdist bdist_wheel

install: clean build
	pip install --upgrade dist/*.whl $(installargs)

install_test_requirements:
	pip install -r requirements-test.txt

install_extra_requirements:
	pip install -r requirements-extra.txt

test: install install_extra_requirements install_test_requirements
	coverage run -m pytest $(pytestopts) $(tests)
	coverage report -m
	coverage xml

test_release_setup:
	twine check dist/*

docs:
	make -C doc html

lint:
	pylint $(package)

reformat:
	black $(package)
	black $(tests)

docker:
	# build
	docker build -f Dockerfile -t $(repo):$(version) .
	# add alternate tags
	docker tag $(repo):$(version) $(repo):latest
	# push to Docker Hub
	docker login && docker push $(repo)

tag:
	git tag $(version)

push_tag:
	git push origin --tags

del_tag:
	git tag -d $(version)

pypi_release:
	twine upload dist/*

release: clean tag
	${MAKE} install test pypi_release push_tag || (${MAKE} del_tag && exit 1)

	# github release
	curl -v -i -X POST \
		-H "Content-Type:application/json" \
		-H "Authorization: token $(token)" \
		https://api.github.com/repos/$(repo)/releases \
		-d '{\
		  "tag_name":"$(version)",\
		  "target_commitish": "master",\
		  "name": "$(version)",\
		  "body": "$(desc)",\
		  "draft": false,\
		  "prerelease": false \
		}'

# build a package with the files needed to run the workflows
workflow:
	mkdir -p dist
	tar -C paper -czf dist/atropos-paper-workflow.tgz \
		workflow/simulated.nf \
		workflow/rnaseq.nf \
		workflow/wgbs.nf \
		workflow/nextflow.config \
		workflow/run-workflows.sh \
		workflow/bin \
		containers/data/simulated/art_profiles.txt \
		containers/tools/tool-names.txt
