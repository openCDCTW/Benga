make-pyc:
	python -m compileall src
	mv $(PWD)/src/__pycache__/*.pyc bin

compile-ui:
	pyside-uic gui/design.ui -o src/view.py

make-venv:
	virtualenv -p python3 venv
	source venv/bin/activate
	pip install -r requirements.txt

help:
	@echo "    clean-pyc"
	@echo "        Remove python artifacts."
	@echo "    clean-build"
	@echo "        Remove build artifacts."
	@echo "    test"
	@echo "        Run py.test"
