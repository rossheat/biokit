import 'dart:io';
import 'package:biokit/constants/strings.dart';
import 'package:meta/meta.dart';

Future<List<Map<String, String>>> readFASTA({String path, String str}) async {
  String contents = path == null ? str : await File(path).readAsString();
  List<String> lines = contents.split('\n');
  int seqCount = 0;

  List<Map<String, String>> fastaMaps = [];
  Map<String, String> currentMap = {};

  for (var line in lines) {
    if (line.startsWith('>')) {
      // starting new line
      if (seqCount != 0) {
        fastaMaps.add(currentMap);
        currentMap = {};
      }
      seqCount++;
      currentMap[kSeq] = '';

      String topLine = line.split('>')[1];
      List<String> topLineList = topLine.split(' ');

      currentMap[kId] = topLineList.first;
      currentMap[kDesc] = topLineList.sublist(1, topLineList.length).join();
    } else {
      currentMap[kSeq] += line;
    }
  }
  fastaMaps.add(currentMap);
  return fastaMaps;
}
