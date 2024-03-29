IMAGE=scl3/task_hii_water


build:
	docker build --no-cache -t $(IMAGE) .

run:
	docker run --rm -it --env-file .env -v `pwd`/src:/app -v `pwd`/.git:/app/.git $(IMAGE) python task.py

shell:
	docker run -it --env-file .env -v `pwd`/src:/app -v `pwd`/.git:/app/.git $(IMAGE) bash

cleanup:
	isort `pwd`/src/*.py
	black `pwd`/src/*.py
	flake8 `pwd`/src/*.py
	mypy `pwd`/src/*.py