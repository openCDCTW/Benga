import os
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase
from rest_framework.test import APIClient
from django.conf import settings


TESTDATA_ROOT = os.path.join(settings.BASE_DIR, 'testdata')


class SequenceTests(APITestCase):
    def setUp(self):
        self.client = APIClient()

    def tearDown(self):
        self.client = None

    def test_create_seq(self):
        """
        Ensure we can create a new Sequence object.
        """
        url = reverse("sequence-list")
        test_file = os.path.join(TESTDATA_ROOT, "vibrio_n16961.fasta")
        with open(test_file, "rb") as file:
            response = self.client.post(url, {"file": file})
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(response.data.keys(), {"id", "created", "file"},
                         "Recieved object does not contain 'id' and 'created' field.")


class TrackedResultsTests(APITestCase):
    def setUp(self):
        self.client = APIClient()
        url = reverse("sequence-list")
        test_file = os.path.join(TESTDATA_ROOT, "vibrio_n16961.fasta")
        with open(test_file, "rb") as file:
            response = self.client.post(url, {"file": file})
        self.seqid = response.data["id"]

    def tearDown(self):
        self.client = None

    def test_create_results(self):
        """
        Ensure we can create a new TrackedResults object.
        """
        url = reverse("results-list")
        test_file = os.path.join(TESTDATA_ROOT, "test.json")
        with open(test_file, "rb") as file:
            data = {"id": self.seqid, "json": file}
            response = self.client.post(url, data)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(response.data.keys(), {"id", "json"},
                         "Recieved object does not contain 'id' and 'json' field.")
        self.assertEqual(str(response.data["id"]), self.seqid, "Inconsistent in 'id' field.")
