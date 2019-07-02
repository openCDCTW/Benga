from rest_framework import mixins, generics
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.http import Http404
from profiling.models import Batch, Sequence, Profile
from profiling.serializers import BatchSerializer, SequenceSerializer, ProfileSerializer
from profiling.tasks import single_profiling, zip_save


class BatchList(generics.ListCreateAPIView):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer


class BatchDetail(mixins.RetrieveModelMixin,
                  mixins.UpdateModelMixin,
                  mixins.DestroyModelMixin,
                  generics.GenericAPIView):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer

    def get(self, request, *args, **kwargs):
        return self.retrieve(request, *args, **kwargs)

    def patch(self, request, *args, **kwargs):
        return self.partial_update(request, *args, **kwargs)

    def delete(self, request, *args, **kwargs):
        return self.destroy(request, *args, **kwargs)


class SequenceList(generics.ListCreateAPIView):
    def get(self, request, format=None):
        sequence = Sequence.objects.all()
        serializer = SequenceSerializer(sequence)
        return Response(serializer.data)

    def post(self, request, format=None):
        serializer = SequenceSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            seq_id = str(serializer.data["id"])
            seq = SequenceDetail.get_object(pk=seq_id)
            batch_id = str(serializer.data["batch_id"])
            (single_profiling.s(seq_id, batch_id,
                                str(seq.database), str(seq.occurrence)) | zip_save.s())()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class SequenceDetail(APIView):
    @classmethod
    def get_object(cls, pk):
        try:
            return Sequence.objects.get(pk=pk)
        except Sequence.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = SequenceSerializer(sequence)
        return Response(serializer.data)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class BatchSequenceDetail(APIView):
    @classmethod
    def get_object(cls, batch_id):
        try:
            return Sequence.objects.get(batch_id=batch_id)
        except Sequence.DoesNotExist:
            raise Http404

    def get(self, request, batch_id, format=None):
        sequence = self.get_object(batch_id)
        serializer = SequenceSerializer(sequence)
        return Response(serializer.data)


class ProfileList(generics.ListCreateAPIView):
    queryset = Profile.objects.all()
    serializer_class = ProfileSerializer


class ProfileDetail(APIView):
    @classmethod
    def get_object(cls, pk):
        try:
            return Profile.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        profile = self.get_object(pk)
        serializer = ProfileSerializer(profile)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        profile = self.get_object(pk)
        serializer = ProfileSerializer(profile, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        profile = self.get_object(pk)
        profile.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)
