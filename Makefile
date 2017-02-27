test = tests
build = python setup.py build_ext -i && \
		python setup.py install && \
		pytest $(test)

install:
	$(call build,)

release:
	# tag
	git tag $(version)
	# build
	$(call build,)
	python setup.py sdist bdist_wheel
	# release
	twine register dist/atropos-$(version).tar.gz
	twine upload dist/atropos-$(version).tar.gz
	git push origin --tags

docs:
	make -C docs api
	make -C docs html

readme:
	pandoc --from=markdown --to=rst --output=README.rst README.md
