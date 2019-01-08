#!/usr/bin/env bash
gunicorn benga.wsgi:application -w 2 -k uvicorn.workers.UvicornWorker