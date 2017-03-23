tests = tests
module = atropos

BUILD = python setup.py build_ext -i && python setup.py install
TEST = py.test $(tests)

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
	twine register dist/$(module)-$(version).tar.gz
	twine upload dist/$(module)-$(version).tar.gz
	git push origin --tags

docs:
	make -C docs api
	make -C docs html

readme:
	pandoc --from=markdown --to=rst --output=README.rst README.md

lint:
	pylint $(module)
