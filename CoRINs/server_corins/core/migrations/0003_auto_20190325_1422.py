# Generated by Django 2.1.7 on 2019-03-25 17:22

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0002_auto_20190325_1413'),
    ]

    operations = [
        migrations.AlterField(
            model_name='arquivo',
            name='files',
            field=models.FileField(upload_to='../up_dj_ring_files/'),
        ),
    ]