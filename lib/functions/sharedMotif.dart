import 'package:biokit/helpers/genCombos.dart';
import 'package:meta/meta.dart';

// Returns the longest common sub-sequence between n sequences.
sharedMotif({@required List<String> seqs}) {
  // Sort to shortest sequence first.
  seqs.sort((a, b) => a.length.compareTo(b.length));

  // Get the shortest sequence.
  String shortestSeq = seqs.first;

  // Generate all possible combinations of shortest sequence.
  List<String> combinations = genCombos(seq: shortestSeq, sorted: true);

  // The longest sub-sequence common to all strings.
  String longComSeq = '';

  // Find the longest combination that is contained in all sequences.
  for (var comb in combinations) {
    bool allMatches = true;
    for (var seq in seqs) {
      if (!seq.contains(comb)) {
        allMatches = false;
        break;
      }
    }
    if (allMatches) {
      longComSeq = comb;
      break;
    }
  }
  return longComSeq;
}
