import os
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase
from rest_framework.test import APIClient
from django.conf import settings


TESTDATA_ROOT = os.path.join(settings.BASE_DIR, 'test', 'data')
