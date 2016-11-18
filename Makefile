build = python setup.py build_ext -i && \
		python setup.py install && \
		nosetests -P tests

install:
	$(call build,)

release:
	# tag
	git tag $(version)
	git push origin --tags
	# build
	$(call build,)
	python setup.py sdist bdist_wheel
	# release
	twine register dist/atropos-$(version).tar.gz
	twine upload dist/atropos-$(version).tar.gz

docs:
	make -C docs api
	make -C docs html

readme:
	pandoc --from=markdown --to=rst --output=README.rst README.md
