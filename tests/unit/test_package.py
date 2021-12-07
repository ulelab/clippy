def test_can_import_package():
    import clip


def test_can_import_interaction_submodule():
    import clip.interaction


def test_package_has_version_string():
    import clip

    assert isinstance(clip.__version__, str)
