#!/usr/bin/env bash
gunicorn benga.asgi:application -w 4 -k uvicorn.workers.UvicornWorker