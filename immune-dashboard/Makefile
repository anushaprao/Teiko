.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt --break-system-packages --quiet

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	python app.py
