def test_can_import_package():
    import clip
    import clip.interaction

def test_package_has_version_string():
    import clip
    assert(isinstance(clip.__version__, str))
