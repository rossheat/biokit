import 'package:meta/meta.dart';

List<String> rfToProtein({@required String aaSeq}) {
  List<String> currentProtein = [];
  List<String> proteins = [];

  for (var aa in aaSeq.split('')) {
    if (aa == 'X') {
      if (currentProtein != []) {
        for (var pro in currentProtein) {
          proteins.add(pro);
        }
        currentProtein = [];
      }
    } else {
      if (aa == 'M') {
        currentProtein.add('');
      }
      for (var i = 0; i < currentProtein.length; i++) {
        currentProtein[i] += aa;
      }
    }
  }
  return proteins;
}
