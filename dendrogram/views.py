from django.http import Http404
from rest_framework import mixins, generics, status
from rest_framework.response import Response
from rest_framework.views import APIView
from django.contrib.sites.models import Site
from django.urls import reverse
from django.core.files import File
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
        filename = request.data["file"].name
        serializer = ProfileSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            batch_id = str(serializer.data["batch_id"])
            linkage = Batch.objects.get(pk=batch_id).prof_num
            prof_num = Batch.objects.get(pk=batch_id).linkage
            url = Site.objects.get_current().domain + reverse("dendrogram-list")
            plot_dendrogram.delay(batch_id, linkage, prof_num,
                                  File(open(serializer.data["file"], "rb")),
                                  filename, url)
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class ProfileDetail(APIView):
    def get_object(self, pk):
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
    def get_object(self, pk):
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
