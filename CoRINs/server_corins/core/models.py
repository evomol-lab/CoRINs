from django.db import models

# Create your models here.
class Arquivo(models.Model):
    files = models.FileField(upload_to='up_dj_ring_files/', blank=False, null=False)