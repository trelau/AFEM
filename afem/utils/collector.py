from collections import defaultdict

__all__ = ["Collector"]


class Collector(defaultdict):
    def __init__(self):
        super(Collector, self).__init__(set)

    @property
    def is_empty(self):
        if self:
            return False
        return True

    @property
    def has_items(self):
        if self:
            return True
        return False

    def add_items(self, key, *items):
        """
        Add items.

        :param key:
        :param items:

        :return:
        """
        for item in items:
            self[key].add(item)

    def get_items(self, key, indx=None):
        """
        Get items.

        :param key:
        :param indx:

        :return:
        """
        if key is None:
            return self.all_items()
        items = list(self[key])
        if indx is None:
            return items
        try:
            return items[indx]
        except IndexError:
            return None

    def all_items(self):
        """
        Get all items.

        :return:
        """
        set_out = set()
        for iset in self.values():
            set_out.update(iset)
        return list(set_out)

    def clear_items(self, key=None):
        """
        Clear items

        :param key:
        :return:
        """
        if key is None:
            self.clear()
            return True
        try:
            self.pop(key)
            return True
        except KeyError:
            return False

    def items_to_key(self):
        """
        Generate a dictionary where the Collector items are the keys and the
        keys are now the values.

        :return:
        """
        items_to_key = dict()
        for key in self.keys():
            for item in self.get_items(key):
                items_to_key[item] = key
        return items_to_key

    def count_items(self, key=None):
        """
        Count items by key.

        :param key:

        :return:
        :rtype: int
        """
        return len(self.get_items(key))

    def remove_items(self, *items):
        """
        Remove items.

        :param items: Items to remove.

        :return:
        """
        status = {}
        for item in items:
            for key in self.keys():
                try:
                    self[key].remove(item)
                    status[item] = True
                except KeyError:
                    status[item] = False
        return status
