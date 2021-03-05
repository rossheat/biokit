import 'package:meta/meta.dart';

dynamic proMotifToRe({@required motifSeq}) {
  String re = '';

  List<String> chars = motifSeq.split('');

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
