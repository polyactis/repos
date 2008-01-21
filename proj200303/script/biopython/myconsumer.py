#!/usr/bin/python2.2

from Bio.Blast import NCBIStandalone
import string
from Bio.ParserSupport import AbstractConsumer

class SequencesExtractor(AbstractConsumer):

    def __init__(self):
        self.sequences_list = []

    def title(self, title_info):
        title_atoms = string.split(title_info)
        new_s = title_atoms[0]

        self.sequences_list.append(new_s)


def extract_sequences(file):
    scanner = NCBIStandalone._Scanner()
    consumer = SequencesExtractor()

    file_to_parse = open(file, 'r')

    
    scanner.feed(file_to_parse, consumer)

    file_to_parse.close()

    return consumer.sequences_list

sequences=extract_sequences("tmp")

print len(sequences)