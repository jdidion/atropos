tests  = tests

BUILD = python setup.py build_ext -i && python setup.py install
TEST  = pytest $(tests)

all:
	$(BUILD)
	$(TEST)

install:
	$(BUILD)

test:
	$(TEST)

release:
	# tag
	git tag $(version)
	# build
	$(BUILD)
	$(TEST)
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
