#!/usr/bin/python
# -*- coding: utf-8 -*-


class BadSettings(Exception):

    """
    One of the following situations may raise this error:
        1. An essential setting wasn't chosen
    """


class OutputCreation(Exception):

    """
    One of the following situations may raise this error:
        1.
        2.

    """


class CatalogueError(Exception):

    """
    One of the following situations may raise this error:
        1. Catalogue couldn't be created
    """


class FolderNotPresent(Exception):
    """Error while execution tried to acces a folder
    One of the following situations may raise this error:
        1. Folder is not present
    """


class FilesNotPresent(Exception):
    """
    One of the following situations may raise this error:
        1. A file is not present
    """


class FolderNotCreated(Exception):
    """
    One of the following situations may raise this error:
        1. Folder couldn't be created
    """


class FullPipelineFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class CleanFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class SplitFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class SextractorFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class ScampFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class FiltFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class RestartFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class ChangeTimeFailed(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class WrongOS(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class AllSameException(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class InvalidScampConfiguration(Exception):
    """
    One of the following situations may raise this error:
        1.
    """


class WrongTicksList(Exception):
    """
    One of the following situations may raise this error:
        1.
    """
