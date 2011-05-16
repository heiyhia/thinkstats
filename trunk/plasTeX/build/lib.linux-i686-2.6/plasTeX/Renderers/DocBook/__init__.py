#!/usr/bin/env python
import re
from plasTeX.Renderers.PageTemplate import Renderer as _Renderer

class DocBook(_Renderer):
    """ Renderer for DocBook documents """
    fileExtension = '.xml'
    imageTypes = ['.png','.jpg','.jpeg','.gif']
    vectorImageTypes = ['.svg']

    def cleanup(self, document, files, postProcess=None):
        res = _Renderer.cleanup(self, document, files, postProcess=postProcess)
        return res

    def processFileContent(self, document, s):
        s = _Renderer.processFileContent(self, document, s)

        # Force XHTML syntax on empty tags
        s = re.compile(r'(<(?:hr|br|img|link|meta|col)\b.*?)\s*/?\s*(>)', 
                       re.I|re.S).sub(r'\1 /\2', s)

        s = re.compile(r'xml:id',re.I).sub(r'id',s)
        #s = re.compile(r'</partintro>\s*<partintro>',re.I).sub(r'',s)

        # replace the first chapter with a preface
        s = re.compile(r'<chapter',re.I).sub(r'<preface', s, count=1)
        s = re.compile(r'</chapter>',re.I).sub(r'</preface>', s, count=1)
        #
        # get rid of duplicate xref gentext
        s = re.compile(r'Chapter.<xref',re.I).sub(r'<xref', s)
        s = re.compile(r'Section.<xref',re.I).sub(r'<xref', s)
        s = re.compile(r'Exercise.<xref',re.I).sub(r'<xref', s)
        s = re.compile(r'Figure.<xref',re.I).sub(r'<xref', s)
        #
        terms = ['figure', 'bookinfo', 'informalfigure', 'programlisting',
                 'informalexample', 'imageobject', 'blockquote', 'sidebar']
        pattern1 = r'<para>\s*(<(%s))' % or_terms(terms)
        pattern2 = r'(</(%s)>)\s*</para>' % or_terms(terms)
        s = re.compile(pattern1, re.I).sub(r'\1',s)
        s = re.compile(pattern2, re.I).sub(r'\1',s)
        #
        s = re.compile(r'<para>([^<]*<title)',re.I).sub(r'\1',s)
        #
        s = re.compile(r'(</mediaobject>)\s*</para>',re.I).sub(r'\1',s)
        #
        s = re.compile(r'<para>\s*(<anchor[^>]*>)\s*</para>',re.I).sub(r'',s)
        #
        s = re.compile(r'(<para>)(\s*<para>)+',re.I).sub(r'\1',s)
        s = re.compile(r'(</para>\s*)+(</para>)',re.I).sub(r'\2',s)
        #
        #pattern1 = r'(<listitem>)\s*<para>'
        #pattern2 = r'</para>\s*(</listitem>)'
        #s = re.compile(pattern1, re.I).sub(r'\1',s)
        #s = re.compile(pattern2, re.I).sub(r'\1',s)
        #
        s = re.compile(r'<para>\s*</para>', re.I).sub(r'', s)
        s = re.compile(r'(<programlisting>)\s*', re.I).sub(r'\1', s)
        #
        return s
    
def or_terms(t):
    return '|'.join(['(%s)' % s for s in t])

Renderer = DocBook
