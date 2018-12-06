#!/usr/bin/env bash
export DATABASE_NAME=""
export DATABASE_USER=""
export DATABASE_PASSWORD=""
gunicorn benga.asgi:application -w 4 -k uvicorn.workers.UvicornWorker