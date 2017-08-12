__all__ = ["order_parts_by_id"]


def order_parts_by_id(parts):
    """
    Order the list of parts by id.

    :param list[afem.structure.entities.Part] parts: The parts.

    :return: List of parts sorted by their ID.
    :rtype: list[afem.structure.entities.Part]
    """
    part_order = [(part.id, part) for part in parts]
    part_order.sort(key=lambda tup: tup[0])
    return [row[1] for row in part_order]
