from rest_framework import mixins, generics
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.http import Http404
from profiling.models import Batch, Sequence, Profile
from profiling.serializers import BatchSerializer, SequenceSerializer,\
    ProfileSerializer, ProfilingSerializer
from profiling.tasks import batch_profiling, profile_and_tree


class BatchList(generics.ListCreateAPIView):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer


class BatchDetail(mixins.RetrieveModelMixin,
                  mixins.DestroyModelMixin,
                  generics.GenericAPIView):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer

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


class Profiling(APIView):
    def post(self, request, format=None):
        serializer = ProfilingSerializer(data=request.data)
        if serializer.is_valid():
            batch_profiling.delay(str(serializer.data["id"]), serializer.data["database"],
                                  serializer.data["occurrence"])
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class ProfilingTree(APIView):
    def post(self, request, format=None):
        serializer = ProfilingSerializer(data=request.data)
        if serializer.is_valid():
            profile_and_tree.delay(str(serializer.data["id"]), serializer.data["database"],
                                   serializer.data["occurrence"])
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
