from django.urls import path

from . import views

urlpatterns = [
    path('profile/', views.ProfileList.as_view(), name="profile-list"),
    path('profile/<uuid:pk>/', views.ProfileDetail.as_view(), name="profile-detail"),
    path('results/', views.TrackedResultsList.as_view(), name="results-list"),
    path('results/<uuid:pk>/', views.TrackedResultsDetail.as_view(), name="results-detail"),
]
