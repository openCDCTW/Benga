from rest_framework import viewsets
from profiling.models import UploadBatch, Sequence, Profile
from profiling.serializers import UploadBatchSerializer, SequenceSerializer, ProfileSerializer

class UploadBatchViewSet(viewsets.ModelViewSet):
    """
        This viewset automatically provides `list` and `detail` actions.
        """
    queryset = UploadBatch.objects.all()
    serializer_class = UploadBatchSerializer


class SequenceViewSet(viewsets.ModelViewSet):
    """
        This viewset automatically provides `list` and `detail` actions.
        """
    queryset = Sequence.objects.all()
    serializer_class = SequenceSerializer


class ProfileViewSet(viewsets.ModelViewSet):
    """
    This viewset automatically provides `list` and `detail` actions.
    """
    queryset = Profile.objects.all()
    serializer_class = ProfileSerializer
