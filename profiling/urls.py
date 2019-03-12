from django.urls import path

from . import views

urlpatterns = [
    path('upload/', views.BatchList.as_view(), name="upload-list"),
    path('upload/<uuid:pk>/', views.BatchDetail.as_view(), name="upload-detail"),
    path('sequence/', views.SequenceList.as_view(), name="sequence-list"),
    path('sequence/<uuid:pk>/', views.SequenceDetail.as_view(), name="sequence-detail"),
    path('profile/', views.ProfileList.as_view(), name="profile-list"),
    path('profile/<uuid:pk>/', views.ProfileDetail.as_view(), name="profile-detail"),
]
