from django.http import Http404
from rest_framework import mixins, generics, status
from rest_framework.response import Response
from rest_framework.views import APIView

from dendrogram.models import Batch, Profile, Dendrogram
from dendrogram.serializers import BatchSerializer, ProfileSerializer, DendrogramSerializer
from dendrogram.tasks import plot_dendrogram


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


class ProfileList(generics.ListCreateAPIView):
    def get(self, request, format=None):
        profile = Profile.objects.all()
        serializer = ProfileSerializer(profile)
        return Response(serializer.data)

    def post(self, request, format=None):
        serializer = ProfileSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            plot_dendrogram.delay(str(serializer.data["batch_id"]))
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


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

    def delete(self, request, pk, format=None):
        profile = self.get_object(pk)
        profile.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class DendrogramList(generics.ListCreateAPIView):
    queryset = Dendrogram.objects.all()
    serializer_class = DendrogramSerializer


class DendrogramDetail(APIView):
    @classmethod
    def get_object(cls, pk):
        try:
            return Dendrogram.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        serializer = DendrogramSerializer(dendrogram)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        serializer = DendrogramSerializer(dendrogram, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        dendrogram = self.get_object(pk)
        dendrogram.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)
