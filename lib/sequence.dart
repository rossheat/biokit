import 'package:biokit/dna.dart';
import 'package:biokit/rna.dart';
import 'package:biokit/strings.dart';
import 'package:biokit/structs.dart';
import 'package:biokit/utils.dart';
import 'errors.dart';

class Sequence {
  late String _seq;
  late String _type;
  late int _len;
  late String name;
  late String id;
  late String desc;

  Sequence(
      {required String seq,
      required String type,
      String name = 'Default name',
      String id = 'Default ID',
      String desc = 'Default description'}) {
    this._type = _validateType(type: type);
    this._seq = _validateSeq(seq: seq);
    this.name = name;
    this.id = id;
    this.desc = desc;
  }

  String _validateType({required String type}) {
    String lType = type.toLowerCase();
    if (![kDNA, kRNA, kPep].contains(lType)) {
      throw (Errors.invalidType(type: lType));
    }
    return lType;
  }

  String _validateSeq({required String seq}) {
    int seqLen = seq.length;
    if (seqLen < 3) {
      throw ('Invalid Sequence Length Error. Sequence has less than three elements.');
    }

    String uSeq = seq.toUpperCase();
    uSeq.split('').asMap().forEach((idx, mon) {
      if (this._type == kDNA) {
        if (!Structs.dnaNucs.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      } else if (this.type == kRNA) {
        if (!Structs.rnaNucs.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      } else {
        if (!Structs.aminoAcids.contains(mon)) {
          throw (Errors.invalidSeq(mon: mon, idx: idx, type: this._type));
        }
      }
    });
    this._len = seqLen;
    return uSeq;
  }

  String operator +(Sequence oSeq) {
    if (this._type != oSeq.type) {
      throw ('Cannot add ${this._type} and ${oSeq.type} sequence.');
    }
    return this._seq + oSeq.seq;
  }

  Map<String, dynamic> info() {
    return {
      'seq': seq,
      'type': type,
      'len': len,
      'name': this.name,
      'id': this.id,
      'desc': this.desc
    };
  }

  @override
  String toString() => info().toString();

  String get seq => this._seq;

  String get type => this._type;

  int get len => this._len;

  /// Calculates the frequency of each monomer.
  Map<String, int> freq() {
    Map<String, int> freqMap = {};
    this._seq.split('').forEach((mon) {
      if (freqMap.containsKey(mon)) {
        freqMap[mon] = freqMap[mon]! + 1;
      } else {
        freqMap[mon] = 1;
      }
    });
    return freqMap;
  }

  // For internal use only as part of the complementary function.
  String reversed({required String seq}) => seq.split('').reversed.join('');

  /// Reverses the sequence.
  String reverse() => this._seq.split('').reversed.join('');

  /// Every possible comnbination of the sequence.
  List<String> combinations({sorted = false}) {
    List<String> listSeq = this._seq.split("");
    List<String> combinations = [];
    for (int i = 0; i < listSeq.length; i++) {
      if (i != listSeq.length - 1) {
        combinations.add(listSeq[i]);
      }
      List<String> temp = [listSeq[i]];
      for (int j = i + 1; j < listSeq.length; j++) {
        temp.add(listSeq[j]);
        combinations.add(temp.join());
      }
    }
    if (sorted) {
      // Sort with longest combination first.
      combinations.sort((b, a) => a.length.compareTo(b.length));
      return combinations;
    }
    return combinations;
  }

  /// Finds the indices (zero-based) of a specified motif.
  Map<String, dynamic> findMotif({required String motif, overlap = true}) {
    List<Map<String, dynamic>> matchData = [];
    Map<String, dynamic> matchMotifMap = {};
    String tempRegexMotif = Utils.motifToRe(motif: motif);
    RegExp regexMotif = overlap ? RegExp('(?=$tempRegexMotif)') : RegExp(tempRegexMotif);
    Iterable<RegExpMatch> allMatches = regexMotif.allMatches(seq);
    for (RegExpMatch match in allMatches) {
      matchData.add({
        kMatch: motif,
        kStartIndex: match.start,
        kEndIndex: match.start + motif.length - 1,
      });
    }
    matchMotifMap[kMatchCount] = allMatches.length;
    matchMotifMap[kMatchIndices] = matchData;
    return matchMotifMap;
  }

  // The number of positional differences.
  int difference({required Sequence oSeq}) {
    if (this._len != oSeq.len) {
      throw ('Sequences must be of the same length to calculate difference.');
    }
    if (this._type != oSeq.type) {
      throw ('Sequences must be of the same type to calculate difference.');
    }

    int differenceCount = 0;
    this._seq.split('').asMap().forEach((idx, mon) {
      if (mon != oSeq.seq[idx]) {
        differenceCount++;
      }
    });
    return differenceCount;
  }

  /// Removes all occurrences of the specified motif.
  String splice({required String motif}) {
    String vMotif = _validateSeq(seq: motif);
    return seq.replaceAll(vMotif, '');
  }

  // The longest shared motif.
  String sharedMotif({required oSeq}) {
    if (this._type != oSeq.type) {
      throw ('Cannot find shared motif between ${this._type} and ${oSeq.type} sequence.');
    }

    // Generate all possible combinations.
    List<String> combos = combinations(sorted: true);

    String longestShared = '';

    // Find the longest combination that is contained in all sequences.
    for (var comb in combos) {
      bool allMatches = true;
      if (!oSeq.seq.contains(comb)) {
        allMatches = false;
        break;
      }
      if (allMatches) {
        longestShared = comb;
        break;
      }
    }
    return longestShared;
  }
}
