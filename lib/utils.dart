import 'dart:convert';
import 'dart:io';

import 'package:biokit/strings.dart';
import 'package:meta/meta.dart';

class Utils {
  static Future<String> uniprotIdToFASTA({@required String uniprotId}) async {
    Uri uri = Uri.parse('http://www.uniprot.org/uniprot/$uniprotId.fasta');

    var request = await HttpClient().getUrl(uri);
    var response = await request.close();

    await for (var contents in response.transform(Utf8Decoder())) {
      return contents;
    }
    return 'Error retrieving protein with uniprot ID $uniprotId';
  }

  static Future<List<Map<String, String>>> readFASTA({String path, String str}) async {
    String contents = path == null ? str : await File(path).readAsString();
    List<String> lines = contents.split('\n');
    int seqCount = 0;

    List<Map<String, String>> fastaMaps = [];
    Map<String, String> currentMap = {};

    for (var line in lines) {
      if (line.startsWith('>')) {
        // Starting new line
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

  static String motifToRe({@required motif}) {
    String re = '';

    List<String> chars = motif.split('');

    bool inBrac = false;
    for (String char in chars) {
      if (char == '[') {
        inBrac = true;
        re += char;
      } else if (char == ']') {
        inBrac = false;
        re += char;
      } else {
        if (inBrac) {
          re += char + '|';
        } else {
          re += char;
        }
      }
    }
    return re.replaceAll('{', '[^').replaceAll('}', ']');
  }
}
