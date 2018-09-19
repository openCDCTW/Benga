from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase
from rest_framework.test import APIClient

class UploadbatchTests(APITestCase):
    def setUp(self):
        self.client = APIClient()

    def tearDown(self):
        self.client = None

    def test_get_batch_id(self):
        """
        Ensure we can create a new UploadBatch object.
        """
        url = reverse("upload-list")
        data = {}
        response = self.client.post(url, data, format='json')
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(response.data.keys(), {"id", "created"},
                         "Recieved object does not contain 'id' and 'created' field.")