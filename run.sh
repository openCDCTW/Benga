#!/usr/bin/env bash
gunicorn benga.asgi:application -w 2 -k uvicorn.workers.UvicornWorker