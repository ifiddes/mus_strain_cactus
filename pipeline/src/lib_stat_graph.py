"""
objects and functions to facilitate reading through the output of the filter
pipeline and deriving statistics from the results.
"""
import math
import os
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml
sys.path.append(
  os.path.join(
    os.path.dirname(os.path.abspath(sys.argv[0])), 'filters'))
import lib_filter


class StatNode(object):
  """ StatNode is used to keep track the counts of various labels as we
  traverse the filter output for a genome.
  """
  def __init__(self):
    self.name = ''
    self.parent = None
    self.heirarchy = []
    self.children = []  # list of StatNode objects
    self.childrenNames = {}  # dict of StatNode objects


class StatCount(object):
  """ StatCount is used to store the counts of occurences for a tag.
  """
  def __init__(self, nodeName):
    self.nodeName = nodeName  # at this level of the heirarchy
    self.nodeTranscripts = 0
    self.nodeTranscriptAnnotations = 0
    self.tagTranscripts = 0
    self.tagTranscriptAnnotations = 0
    self.children = []  # more StatCounts


def newElement(parent, tag):
  """
  """
  if tag.startswith('^'):
    tag = tag[1:] + '_prexisting'
  e = ET.SubElement(parent, tag)
  e.attrib['transcripts'] = '0'
  e.attrib['transcript_annotations'] = '0'
  return e


def newStatGraph():
  """ Create and return the initial stat graph object.
  """
  e = ET.Element('stats')
  n = newElement(e, 'ok')
  n = newElement(e, 'not_ok')
  return ET.ElementTree(e)


def buildStatGraph(args):
  """ Read the transcript bed files and return a StatNode based graph
  """
  g = newStatGraph()
  r = g.getroot()
  transcripts = lib_filter.getTranscripts(args.geneCheck, args.geneCheckDetails)
  for t in transcripts:
    if isOk(t.annotations):
      recordOk(r, t)
    else:
      recordNotOk(r, t)
  return g


def recordStatGraph(g, path):
  """ record a stat graph G to a file located in PATH.
  """
  # g = ET.ElementTree(g)
  g.write(path, xml_declaration=True, encording='utf-8', method='xml')


def readStatGraph(xml):
  """ read a stat graph located at path PATH.
  """
  return ET.parse(xml)


def recordOk(root, transcript):
  """ Record in the graph that this transcript is ok.
  """
  ok = root.find('ok')
  addOne(ok, 'transcripts')
  # flatten the annotations so that there are no duplicates
  flattenedTree = flattenAnnotations(transcript.annotations)
  flatTree = flattenedTree
  depthFirstAddOne(flatTree, ok, isRoot=True)

  # walk the full annotation set as it appears
  for ta in transcript.annotations:
    addOne(ok, 'transcript_annotations')
    prev = ok
    for label in ta.labels:
      label = cleanLabel(label)
      e = prev.find(label)
      if e is None:
        raise RuntimeError('Unanticipated tag discovered %s:%s'
                           % (prev.tag, label))
      addOne(e, 'transcript_annotations')
      prev = e


def addOne(tag, at):
  """ For TAG, treat the attribute AT as an int and add one to it.
  """
  tag.attrib[at] = str(int(tag.attrib[at]) + 1)


def cleanLabel(label):
  """ clean up any weird text that would break xml name tag formatting.
  """
  if label.startswith('^'):
    label = label[1:] + '_prexisting'
  return label


def flattenAnnotations(transcriptAnnotations):
  """ for a list of TranscriptAnnotations, create a list that contains
  no duplicate heirarchies.
  """
  root = StatNode()
  root.name = 'not_ok'
  for ta in transcriptAnnotations:
    pos = root
    # descend the implicit heirarchy of the labels and build a tree
    for label in ta.labels:
      label = cleanLabel(label)
      if label not in pos.childrenNames:
        n = StatNode()
        n.name = label
        n.parent = pos
        n.heirarchy = pos.heirarchy[:]
        n.heirarchy.append(n.name)
        pos.childrenNames[label] = n
        pos.children.append(n)
      pos = pos.childrenNames[label]
  return root


def depthFirstAddOne(node, notOkTree, isRoot=False):
  """ descend the tree in a depth-first pre-order fashion and record in the
  notOkTree one 'transcript' per node.
  """
  if not isRoot:
    e = notOkTree
    prev = e
    for label in node.heirarchy:
      label = cleanLabel(label)
      e = prev.find(label)
      if e is None:
        if prev.tag == 'not_ok' and label.startswith('hasBadCopies'):
          print node.name, node.parent.name, map(lambda x: x.name, node.children)
          raise RuntimeError('dfao adding %s:%s, %s' %
                             (prev.tag, label, str(node.heirarchy)))
        e = newElement(prev, label)
      prev = e
    addOne(e, 'transcripts')
  for c in node.children:
    depthFirstAddOne(c, notOkTree)


def recordNotOk(root, transcript):
  """ For the TranscriptAnnotations in the transcript, make a
  record in the graph
  """
  no = root.find('not_ok')
  addOne(no, 'transcripts')

  # flatten the annotations so that there are no duplicates
  flattenedTree = flattenAnnotations(transcript.annotations)
  flatTree = flattenedTree
  depthFirstAddOne(flatTree, no, isRoot=True)

  # walk the full annotation set as it appears
  for ta in transcript.annotations:
    addOne(no, 'transcript_annotations')
    prev = no
    for label in ta.labels:
      label = cleanLabel(label)
      e = prev.find(label)
      if e is None:
        raise RuntimeError('Unanticipated tag discovered %s:%s'
                           % (prev.tag, label))
      addOne(e, 'transcript_annotations')
      prev = e


def isOk(annots):
  """ given a list of TranscriptAnnotations, return True if data is OK.
  """
  if annots == []:
    return True
  for a in annots:
    for label in a.labels:
      if (label != 'hasOkCopies' and
          label != 'hasBadCopies' and
          not label.startswith('count_')):
        return False
  return True


def reportCounts(counts):
  """ Given a dict of counts (key: category, value: integer count, print it out.
  """
  max_str = 0
  max_digit = 0
  for c in counts:
    if max_str <= len(c):
      max_str = len(c) + 1
    try:
      v = round(math.log10(counts[c]))
    except ValueError:
      v = 1
    if max_digit < v:
      max_digit = int(v)
  print '%*s %*d (%.3f)' % (max_str, 'total', max_digit, counts['total'], 1.0)
  print '%*s %*d (%.3f)' % (max_str, 'ok', max_digit, counts['ok'],
                            float(counts['ok']) / counts['total'])
  print '%*s %*d (%.3f)' % (max_str, 'not ok', max_digit, counts['not ok'],
                            float(counts['not ok']) / counts['total'])
  order = sorted(counts, key=lambda x: counts[x], reverse=True)
  for cat in order:
    if cat in ['total', 'ok', 'not ok']:
      continue
    print '%*s %*d (%.3f)' % (max_str, cat, max_digit, counts[cat],
                               float(counts[cat]) / counts['total'])


def graphToStatCount(node, tag):
  """ Given a stat graph NODE, return a StatCount of similar structure.
  """
  s = StatCount(node.tag)
  if node.tag != 'stats':
    s.nodeTranscripts = int(node.attrib['transcripts'])
    s.nodeTranscriptAnnotations = int(node.attrib['transcript_annotations'])
    if node.tag == tag:
      s.tagTranscripts = s.nodeTranscripts
      s.tagTranscriptAnnotations = s.nodeTranscriptAnnotations
  for child in node:
    s.children.append(graphToStatCount(child, tag))
  if node.tag == 'stats':
    for child in s.children:
      s.nodeTranscripts += child.nodeTranscripts
      s.nodeTranscriptAnnotations += child.nodeTranscriptAnnotations
  return s


def statCountContainsTag(node, tag):
  """ Walk the entire NODE tree and if any tags match TAG, return True.
  """
  if node.nodeName == tag:
    return True
  if node.children == []:
    return False
  c = False
  for child in node.children:
    if statCountContainsTag(child, tag):
      return True
  return False


def pruneStatCountBranches(node, tag):
  """ Given label TAG and StatCount NODE, remove branches without TAG leafs.
  """
  if node.nodeName == tag:
    node.children = []  # we don't care about anything beyond the tag.
  newChildren = []
  for child in node.children:
    if statCountContainsTag(child, tag):
      newChildren.append(child)
  node.children = newChildren
  for child in node.children:
    pruneStatCountBranches(child, tag)


def sendUpStatCountTagCounts(node, tag):
  """ Given label TAG and StatCount NODE, update counts for nodes leaf to root.
  """
  def pushUp(node):
    t = 0
    ta = 0
    for child in node.children:
      tc, tac = pushUp(child)
      t += tc
      ta += tac
    node.tagTranscripts += t
    node.tagTranscriptAnnotations += ta
    return node.tagTranscripts, node.tagTranscriptAnnotations
  pushUp(node)


def getTagStats(graph, tag):
  """ Given a stat graph GRAPH and a tag name TAG, return stats from GRAPH.
  """
  r = graph.getroot()
  s = graphToStatCount(r, tag)
  pruneStatCountBranches(s, tag)
  sendUpStatCountTagCounts(s, tag)
  return s


def reportTagStats(stats, lower):
  """ Given some tag statistics STATS, report the data.
  """
  level = 0
  increment = 2
  def printTree(level, t, lower):
    if t.tagTranscripts > lower:
      print('%s%s %10d, %5d, %5d (%.2f%%), %5d (%.2f%%)' %
            (' ' * level, t.nodeName,
             t.nodeTranscripts, t.nodeTranscriptAnnotations,
             t.tagTranscripts, 100. * t.tagTranscripts / t.nodeTranscripts,
             t.tagTranscriptAnnotations,
             100. * t.tagTranscriptAnnotations / t.nodeTranscriptAnnotations))
    for c in t.children:
      printTree(level + increment, c, lower)
  printTree(level, stats, lower)
