import 'package:biokit/sequence.dart';
import 'package:biokit/strings.dart';
import 'package:biokit/structs.dart';
import 'package:meta/meta.dart';

class Nucleotides extends Sequence {
  Nucleotides({@required String seq, @required String type}) : super(seq: seq, type: type);

  /// The frequency of each nucleotide.
  Map<String, int> freq() => super.freq();

  /// Translates the nucleotides to amino acids.
  Map<String, dynamic> translate({revComp = false, startIdx = 0}) {
    String seq = revComp ? complementary(rev: true) : this.seq;

    String aaSeq = '';
    for (var i = startIdx; i < seq.length - 2; i += 3) {
      String codon = seq.substring(i, i + 3);
      aaSeq += this.type == kDNA ? Structs.dnaCodonToAA[codon] : Structs.rnaCodonToAA[codon];
    }
    return {kAASeq: aaSeq, 'nucCount': seq.length - startIdx - 1, 'aaCount': aaSeq.length};
  }

  /// Calculates the frequency of each codon for a specified amino acid.
  /// Searches in batches of three nucleotides.
  /// If a codon is present but does not appear in a batch, it will not be detected.
  /// For example in ATGTCATGC, only 1 ATG is detected.
  Map<String, int> codonFreq({@required String aminoAcid}) {
    Map<String, int> codonFreqMap = {};
    for (var i = 0; i < this.len - 2; i += 3) {
      String codon = this.seq.substring(i, i + 3);
      String fetchedAminoAcid =
          this.type == kDNA ? Structs.dnaCodonToAA[codon] : Structs.rnaCodonToAA[codon];
      if (fetchedAminoAcid == aminoAcid.toUpperCase()) {
        codonFreqMap[codon] == null ? codonFreqMap[codon] = 1 : codonFreqMap[codon]++;
      }
    }
    return codonFreqMap;
  }

  /// The complementary strand.
  String complementary({bool rev = false}) {
    String compSeq = this
        .seq
        .split('')
        .map((nuc) => this.type == kDNA ? Structs.dnaCompNucs[nuc] : Structs.rnaCompNucs[nuc])
        .join();
    return rev ? super.reversed(seq: compSeq) : compSeq;
  }

  /// The percentage of nucleotides which are either Guanine or Cytosine.
  double gcContent() {
    int gcCount = 0;
    this.seq.split('').forEach((nuc) {
      if (nuc == 'G' || nuc == 'C') {
        gcCount++;
      }
    });
    return num.parse(
      ((gcCount / this.len) * 100).toStringAsFixed(2),
    );
  }

  /// Generates 6 ORFs, 3 forward, 3 reversed.
  List<String> readingFrames() {
    List<String> readingFrames = [];

    for (var i = 0; i < 3; i++) {
      readingFrames.add(translate(startIdx: i)[kAASeq]);
      readingFrames.add(translate(revComp: true, startIdx: i)[kAASeq]);
    }
    return readingFrames;
  }

  List<String> readingFrameToProteins({@required String aaSeq}) {
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

  List<String> proteins({bool unique = false}) {
    List<String> frames = readingFrames();
    List<String> proteins = [];
    for (var frame in frames) {
      List<String> tempProteins = readingFrameToProteins(aaSeq: frame);
      if (proteins != []) {
        proteins.addAll(tempProteins);
      }
    }
    proteins.sort((b, a) => a.length.compareTo(b.length));
    if (unique) {
      return proteins.toSet().toList();
    }
    return proteins;
  }
}
