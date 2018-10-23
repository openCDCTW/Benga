from rest_framework import mixins, generics
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.http import Http404
from profiling.models import UploadBatch, Sequence, Profile, Dendrogram
from profiling.serializers import UploadBatchSerializer, SequenceSerializer,\
    ProfileSerializer, DendrogramSerializer
from profiling.tasks import do_profiling


class UploadBatchList(generics.ListCreateAPIView):
    queryset = UploadBatch.objects.all()
    serializer_class = UploadBatchSerializer


class UploadBatchDetail(mixins.RetrieveModelMixin,
                        mixins.DestroyModelMixin,
                        generics.GenericAPIView):
    queryset = UploadBatch.objects.all()
    serializer_class = UploadBatchSerializer

    def get(self, request, *args, **kwargs):
        return self.retrieve(request, *args, **kwargs)

    def delete(self, request, *args, **kwargs):
        return self.destroy(request, *args, **kwargs)


class SequenceList(generics.ListCreateAPIView):
    queryset = Sequence.objects.all()
    serializer_class = SequenceSerializer


class SequenceDetail(APIView):
    def get_object(self, pk):
        try:
            return Sequence.objects.get(pk=pk)
        except Sequence.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = SequenceSerializer(sequence)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = SequenceSerializer(sequence, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class ProfileList(generics.ListCreateAPIView):
    queryset = Profile.objects.all()
    serializer_class = ProfileSerializer


class ProfileDetail(APIView):
    def get_object(self, pk):
        try:
            return Profile.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = ProfileSerializer(sequence)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = ProfileSerializer(sequence, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class DendrogramList(generics.ListCreateAPIView):
    queryset = Dendrogram.objects.all()
    serializer_class = DendrogramSerializer


class DendrogramDetail(APIView):
    def get_object(self, pk):
        try:
            return Dendrogram.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = DendrogramSerializer(sequence)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = DendrogramSerializer(sequence, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class Profiling(APIView):
    def get_object(self, batch_id):
        try:
            return UploadBatch.objects.get(pk=batch_id)
        except UploadBatch.DoesNotExist:
            return Http404

    def post(self, request, batch_id, format=None):
        batch = self.get_object(batch_id)
        profile = Profile(id=batch)
        serializer = ProfileSerializer(profile, data=request.data)
        if serializer.is_valid():
            do_profiling.delay(batch_id, serializer.data["database"], serializer.data["occr_level"])
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
