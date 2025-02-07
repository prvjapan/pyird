
def update_attr(instance, **kwargs):
    """ Update the attributes using kwargs

    Args:
        instance: Python Class 
        kwargs: keys and values to be updated 
    """
    for key, value in kwargs.items():
        if hasattr(instance, key):
            setattr(instance, key, value)
        else:
            raise AttributeError(f"Class {instance.__class__.__name__} has no attribute named '{key}'")
