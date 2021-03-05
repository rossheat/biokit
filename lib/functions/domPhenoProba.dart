import 'package:meta/meta.dart';

double domPhenoProba(
    {@required int countHomoDom,
    @required int countHet,
    @required int countHomoRec}) {
  int totalOrganisms = countHomoDom + countHet + countHomoRec;

  // Heterozygous * Heterozygous
  double HETHET =
      ((countHet / totalOrganisms) * ((countHet - 1) / (totalOrganisms - 1))) *
          0.25;
  // Heterozygous * Homozygous Recessive
  double HETHR =
      ((countHet / totalOrganisms) * (countHomoRec / (totalOrganisms - 1)));

  // Homozygous Recessive * Homozygous Recessive
  double HRHR = (countHomoRec / totalOrganisms) *
      ((countHomoRec - 1) / (totalOrganisms - 1));

  // 1 - (probability of a recessive phenotype)
  double dominantPhenotypeProbability = 1 - (HETHET + HETHR + HRHR);
  return dominantPhenotypeProbability;
}
