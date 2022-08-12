using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace MandatySim.Con
{
    class Program
    {
        private static readonly RandomNumberGenerator rng = new();

        private const int MandatesTotal = 200;
        private const int IterationCount = 30000;

        private static void Main(string[] args)
        {
            Dictionary<Region, int> citizens = new()
            {
                { new("PR"), 1070032 },
                { new("JC"), 511000 },
                { new("JM"), 949258 },
                { new("KV"), 235018 },
                { new("VY"), 405163 },
                { new("HK"), 438969 },
                { new("LI"), 348860 },
                { new("MS"), 952974 },
                { new("OL"), 501969 },
                { new("PA"), 414895 },
                { new("PL"), 472331 },
                { new("SČ"), 1081557 },
                { new("ÚS"), 644856 },
                { new("ZL"), 464704 },
            };

            var expectedTurnoutGlobal = new RandomVariableCharacteristics(0.595, 5);
            var globalResultsStdDevBig = 6;
            var globalResultsStdDevSmall = 9;
            var expectedResultsGlobal = new Dictionary<Party, RandomVariableCharacteristics>
            {
                { new("Pir+STAN"), new(0.269, globalResultsStdDevBig) },
                { new("ANO"), new(0.213, globalResultsStdDevBig) },
                { new("SPOLU"), new(0.200, globalResultsStdDevBig) },
                { new("SPD"), new(0.112, globalResultsStdDevBig) },
                { new("KSČM"), new(0.055, globalResultsStdDevSmall) },
                { new("ČSSD"), new(0.042, globalResultsStdDevSmall) },
                { new("Přísaha"), new(0.034, globalResultsStdDevSmall) },
                { new("Trik"), new(0.032, globalResultsStdDevSmall) },
                { new("Zelení"), new(0.020, globalResultsStdDevSmall) },
                { new("Jiné"), new(0.037, globalResultsStdDevSmall) },
            };
            Dictionary<Region, RandomVariableCharacteristics> expectedTurnout = citizens.ToDictionary(
                ckv => ckv.Key,
                ckv => expectedTurnoutGlobal
            );
            Dictionary<Region, Dictionary<Party, RandomVariableCharacteristics>> expectedResults = citizens.ToDictionary(
                ckv => ckv.Key,
                ckv => expectedResultsGlobal
            );
            var partyLimitType = expectedResultsGlobal.ToDictionary(
                pkv => pkv.Key,
                pkv => pkv.Key.ShortName switch
                {
                    "Pir+STAN" => 2,
                    "SPOLU" => 3,
                    "Trik" => 3,
                    _ => 1
                }
            );
            var possibleCoalitions = new Dictionary<String, HashSet<Party>>
            {
                { "ANO+Piráti+STAN", new HashSet<Party> { new("ANO"), new("Pir+STAN") } },
                { "Piráti+STAN+SPOLU", new HashSet<Party> { new("Pir+STAN"), new("SPOLU") } },
                { "Piráti+STAN+SPOLU+ČSSD", new HashSet<Party> { new("Pir+STAN"), new("SPOLU"), new("ČSSD") } },
                { "ANO+SPOLU", new HashSet<Party> { new("ANO"), new("SPOLU") } },
                { "ANO+SPOLU+ČSSD", new HashSet<Party> { new("ANO"), new("SPOLU"), new("ČSSD") } },
                { "Piráti+STAN+SPOLU+Přísaha", new HashSet<Party> { new("Pir+STAN"), new("SPOLU"), new("Přísaha") } },
                { "ANO+SPOLU+Přísaha", new HashSet<Party> { new("ANO"), new("SPOLU"), new("Přísaha") } },
                { "ANO+SPD+KSČM", new HashSet<Party> { new("ANO"), new("SPD"), new("KSČM") } },
                { "ANO+SPD+KSČM+ČSSD", new HashSet<Party> { new("ANO"), new("SPD"), new("KSČM"), new("ČSSD") } },
                { "ANO+SPD+KSČM+Přísaha", new HashSet<Party> { new("ANO"), new("SPD"), new("KSČM"), new("Přísaha") } },
                { "ANO+SPD+Přísaha", new HashSet<Party> { new("ANO"), new("SPD"), new("Přísaha") } },
                { "ANO+SPD", new HashSet<Party> { new("ANO"), new("SPD") } },
            };

            var sums = expectedResultsGlobal.ToDictionary(erg => erg.Key, erg => new AtomicInteger(0));
            var minimums = expectedResultsGlobal.ToDictionary(erg => erg.Key, erg => new AtomicInteger(Int32.MaxValue));
            var maximums = expectedResultsGlobal.ToDictionary(erg => erg.Key, erg => new AtomicInteger(Int32.MinValue));
            var partyInParliament = expectedResultsGlobal.ToDictionary(erg => erg.Key, erg => new AtomicInteger(0));

            var coalitionsWorkingCount = possibleCoalitions.ToDictionary(pckv => pckv.Key, pckv => new AtomicInteger(0));
            var coalitionsAllRequiredCount = possibleCoalitions.ToDictionary(pckv => pckv.Key, pckv => new AtomicInteger(0));

            Console.Write("Simulating...");
            RunInParallel(IterationCount, iter =>
            {
                if (iter % 100 == 0) Console.Write('.');

                // 1. generate VotingResults
                Dictionary<Region, Dictionary<Party, int>> results = GenerateRandomResults(citizens, expectedTurnout, expectedResults);

                // 2. compute mandates
                Dictionary<Party, int> mandates = Evaluate(results, partyLimitType);

                foreach (var party in expectedResultsGlobal.Keys)
                {
                    mandates.TryGetValue(party, out var count);
                    sums[party].Add(count);
                    minimums[party].MinWith(count);
                    maximums[party].MaxWith(count);
                    if (count > 0) partyInParliament[party].Increment();
                }

                foreach (var coalition in possibleCoalitions)
                {
                    var count = mandates.Where(mkv => coalition.Value.Contains(mkv.Key)).Sum(mkv => mkv.Value);
                    if (count > MandatesTotal / 2)
                    {
                        coalitionsWorkingCount[coalition.Key].Increment();
                        if (coalition.Value.All(p => mandates.TryGetValue(p, out var pm) && pm > 0)) coalitionsAllRequiredCount[coalition.Key].Increment();
                    }
                }
            });
            Console.WriteLine();

            foreach (var pkv in sums.OrderByDescending(m => m.Value.Value).Where(p => maximums[p.Key].Value > 0))
            {
                Console.WriteLine($"{pkv.Key.ShortName}\t{pkv.Value.Value / (double) IterationCount:N2}\t({minimums[pkv.Key]}–{maximums[pkv.Key]}: {partyInParliament[pkv.Key].Value * 100.0 / IterationCount:N2} %)");
            }

            Console.WriteLine();

            foreach (var ckv in coalitionsWorkingCount.OrderByDescending(c => coalitionsAllRequiredCount[c.Key].Value))
            {
                Console.WriteLine($"{ckv.Key}\t{ckv.Value.Value * 100.0 / IterationCount:N2} % ({coalitionsAllRequiredCount[ckv.Key].Value * 100.0 / IterationCount:N2} %)");
            }
        }

        private static void ShowRandomHistogram()
        {
            var buckets = new int[80];
            for (int i = 0; i < 10000; ++i)
            {
                // var n = (int) Math.Round(MathNet.Numerics.Distributions.LogNormal.Sample(rng,  Math.Log(40), 1));
                // var n = (int) Math.Round(GenerateLogNormalRandomNumber(new(Math.Log(40), 1)));
                var n = GenerateBinomialRandomNumber(8000, 0.4, 1) / 100;
                if (n >= 0 && n < buckets.Length) ++buckets[n];
            }
            var cdf = Enumerable.Range(0, 80).Select(k => GenerateBinomialCdf(80, 0.4, 1 * k)).ToArray();
            var max = buckets.Max();
            for (int row = 20; row >= 0; --row)
            {
                for (int i = 0; i < buckets.Length; ++i)
                {
                    var b = buckets[i];
                    var c = cdf[i];
                    var showB = 21 * b / max >= row;
                    var showC = 21 * c >= row;
                    Console.Write(showB ? (showC ? '#' : '*') : (showC ? '+' : '.'));
                }
                Console.WriteLine();
            }
        }

        private static void RunInParallel(int iterationCount, Action<int> action)
        {
#if DEBUG
            for (var iter = 0; iter < iterationCount; ++iter)
            {
                action(iter);
            }
#else
            Parallel.For(0, iterationCount, action);
#endif
        }

        private static Dictionary<Party, int> Evaluate(Dictionary<Region, Dictionary<Party, int>> results, Dictionary<Party, int> partyLimitType)
        {
            var regionVotes = results.ToDictionary(rkv => rkv.Key, rkv => rkv.Value.Sum(pkv => pkv.Value));
            var totalVotes = regionVotes.Values.Sum();

            // § 48 Určení počtu poslanců volených ve volebních krajích
            //
            // (1) Na základě výsledků hlasování převzatých z volebních okrsků a zvláštních volebních okrsků u pověřených obecních úřadů podle § 43 zjistí Český statistický úřad celkový počet platných hlasů,
            // které byly ve všech volebních krajích odevzdány pro všechny politické strany, politická hnutí a koalice, a vydělí ho počtem volených poslanců podle § 24. Číslo takto vypočtené a zaokrouhlené
            // na jednotky je republikovým mandátovým číslem.
            var republicMandateNumber = (totalVotes + 100) / MandatesTotal;

            // (2) Republikovým mandátovým číslem se dělí celkový počet platných hlasů odevzdaných v každém volebním kraji.
            var regionMandateDivision = regionVotes.ToDictionary(rkv => rkv.Key, rkv => (quot: Math.DivRem(rkv.Value, republicMandateNumber, out int rem), rem));
            // Celé číslo takto vypočtené je počtem mandátů, které připadají jednotlivým volebním krajům.
            var regionMandates = regionMandateDivision.ToDictionary(rkv => rkv.Key, rkv => rkv.Value.quot);
            var remaining = MandatesTotal - regionMandates.Values.Sum();
            if (remaining > 0)
            {
                // (3) Nebyly-li takto rozděleny všechny mandáty, připadnou zbylé mandáty postupně volebním krajům, které vykazují největší zbytek dělení. 
                foreach (var remainderGroup in regionMandateDivision.GroupBy(rkv => rkv.Value.rem).OrderByDescending(g => g.Key))
                {
                    var group = remainderGroup.ToList();
                    while (remaining > 0 && group.Count > 0)
                    {
                        var count = group.Count;
                        // Při rovnosti zbytků rozhoduje los.
                        var chosen = rng.Next(count);
                        var (chosenRegion, _) = group[chosen];
                        group.RemoveAt(chosen);
                        regionMandates[chosenRegion]++;
                        --remaining;
                    }
                    if (remaining == 0) break;
                }
            }

            // § 49 Postup politických stran, politických hnutí a koalic do prvního skrutinia
            //
            // (1) Český statistický úřad zjistí,
            //     a) které politické strany nebo politická hnutí získaly méně než 5 procent z celkového počtu platných hlasů,
            //     b) které koalice složené ze 2 politických stran, popřípadě politických hnutí získaly méně než 8 procent z celkového počtu platných hlasů,
            //     c) které koalice složené ze 3 a více politických stran, popřípadě politických hnutí získaly méně než 11 procent z celkového počtu platných hlasů.
            //
            // (4) Pokud do skrutinia nepostupují politické strany, politická hnutí nebo koalice v počtu podle odstavce 3, Český statistický úřad sníží u
            //     a) politických stran nebo politických hnutí hranici 5 procent na hranici 4 procent z celkového počtu platných hlasů,
            //     b) koalic podle odstavce 1 písm. b) hranici 8 procent na hranici 7 procent z celkového počtu platných hlasů,
            //     c) koalic podle odstavce 1 písm. c) hranici 11 procent na hranici 10 procent z celkového počtu platných hlasů.
            // (5) Nebude-li ani postupem podle odstavce 4 dosaženo postupu do skrutinia v počtu podle odstavce 3, sníží Český statistický úřad hranici o další procento.
            Dictionary<Party, int> partyTotalVotes = null;
            for (var limitIteration = 0; limitIteration < 12; ++limitIteration)
            {
                var limits = Enumerable.Range(0, 3).Select(type => (totalVotes * Math.Max(5 + type * 3 - limitIteration, 0) + 50) / 100).ToArray();
                partyTotalVotes = results.Values.SelectMany(parties => parties).GroupBy(kv => kv.Key).ToDictionary(g => g.Key, g => g.Sum(kv => kv.Value));

                var failingParties = new HashSet<Party>(partyTotalVotes.Where(pkv => pkv.Value < limits[partyLimitType[pkv.Key] - 1]).Select(pkv => pkv.Key));
                // „Jiné“ nikdy nebudou mít dost hlasů
                failingParties.Add(new Party("Jiné"));

                // (3) Český statistický úřad zjistí, zda do skrutinia postupují alespoň
                //     a) 2 koalice,
                //     b) 1 koalice a 1 politická strana nebo politické hnutí, nebo
                //     c) 2 politické strany nebo politická hnutí.                
                if (partyTotalVotes.Keys.Count() - failingParties.Count >= 2)
                {
                    // (2) Při dalším zjišťování volebních výsledků a přidělování mandátů se již k politickým stranám, politickým hnutím a koalicím podle odstavce 1 a hlasům pro ně odevzdaným nepřihlíží.
                    foreach (var failingParty in failingParties)
                    {
                        partyTotalVotes[failingParty] = 0;
                        foreach (var regionResults in results.Values)
                        {
                            regionResults[failingParty] = 0;
                        }
                    }
                    break;
                }
            }

            var partyMandates = new Dictionary<Party, int>();
            var totalAssignedMandates = 0;
            // Do druhého skrutinia se převádějí zbytky hlasů jednotlivých politických stran, politických hnutí a koalic pro ně odevzdaných, a nedostala-li politická
            // strana, politické hnutí nebo koalice v prvním skrutiniu ani jeden mandát, potom všechny hlasy pro ni odevzdané.
            var partyRemainders = new Dictionary<Party, int>();

            // § 50 První skrutinium
            // (1) V prvním skrutiniu se rozdělují mandáty v rámci volebních krajů.
            foreach (var rkv in results)
            {
                // (2) Součet platných hlasů odevzdaných ve volebním kraji pro politické strany, politická hnutí a koalice, které postoupily do prvního skrutinia, se vydělí počtem mandátů, které byly tomuto
                // volebnímu kraji přiděleny, zvětšeným o dvě; číslo takto vypočtené a zaokrouhlené na jednotky je krajským volebním číslem.
                var regionMandateCount = regionMandates[rkv.Key];
                var fixedMandateCount = regionMandateCount + 2;
                var regionMandateNumber = (rkv.Value.Values.Sum() + fixedMandateCount / 2) / fixedMandateCount;

                // (3) Celkový počet platných hlasů, který obdržela politická strana, politické hnutí nebo koalice v rámci volebního kraje, se dělí krajským volebním číslem a politické straně, politickému
                // hnutí nebo koalici se přikáže tolik mandátů, kolikrát je krajské volební číslo obsaženo v celkovém počtu platných hlasů, které tato politická strana, politické hnutí nebo koalice získala.
                var partyMandateDivision = rkv.Value.Where(pkv => pkv.Value > 0).ToDictionary(pkv => pkv.Key, pkv => (quot: Math.DivRem(pkv.Value, regionMandateNumber, out var rem), rem));
                var assignedMandates = 0;
                foreach (var pkv in partyMandateDivision)
                {
                    partyMandates.TryGetValue(pkv.Key, out var prevMandates);
                    partyMandates[pkv.Key] = prevMandates + pkv.Value.quot;
                    partyRemainders.TryGetValue(pkv.Key, out var prevRemainder);
                    partyRemainders[pkv.Key] = prevRemainder + pkv.Value.rem;
                    assignedMandates += pkv.Value.quot;
                    totalAssignedMandates += pkv.Value.quot;
                }

                // (4) Bylo-li takto rozděleno více mandátů, než se mělo přidělit podle § 48 odst. 2, odečtou se přebývající mandáty postupně těm politickým stranám, politickým hnutím nebo koalicím, které
                // ve volebním kraji vykázaly nejmenší zbytek dělení.
                while (assignedMandates > regionMandateCount)
                {
                    // [přidáno: nebudeme odebírat těm, kterým jsme nedali ani jeden mandát!]
                    foreach (var unluckyPartyRemGroup in partyMandateDivision.Where(pmd => partyMandates[pmd.Key] > 0).GroupBy(pmd => pmd.Value.rem).OrderBy(g => g.Key))
                    {
                        // Při stejném zbytku dělení se mandát odečte politické straně, politickému hnutí nebo koalici, která získala ve volebním kraji menší počet hlasů;
                        foreach (var unluckyPartyVoteGroup in unluckyPartyRemGroup.GroupBy(g => rkv.Value[g.Key]).OrderBy(g => g.Key))
                        {
                            var unluckyGroup = unluckyPartyVoteGroup.ToList();
                            // je-li i tak počet platných hlasů stejný, rozhodne los.
                            var chosen = rng.Next(unluckyGroup.Count);
                            var unluckyParty = unluckyGroup[chosen].Key;
                            --partyMandates[unluckyParty];
                            Debug.Assert(partyMandates[unluckyParty] >= 0);
                            --assignedMandates;
                            --totalAssignedMandates;
                            partyRemainders[unluckyParty] += regionMandateNumber;

                            if (assignedMandates == regionMandateCount) break;
                        }

                        if (assignedMandates == regionMandateCount) break;
                    }
                }
                Debug.Assert(assignedMandates <= regionMandateCount);
            }

            // § 51 Druhé skrutinium
            // (1) Všechny mandáty, které nebyly přikázány v prvním skrutiniu, se přikáží ve druhém skrutiniu.
            var remainingMandates = MandatesTotal - totalAssignedMandates;
            if (remainingMandates > 0)
            {
                // (2) Ve druhém skrutiniu se sečtou zbytky hlasů jednotlivých politických stran, politických hnutí a koalic. Tento součet se vydělí počtem mandátů, které nebyly v prvním skrutiniu přikázány,
                // zvětšeným o jednu. Číslo takto vypočtené a zaokrouhlené na jednotky je republikovým volebním číslem.
                var totalRemainders = partyRemainders.Values.Sum();
                var fixedRemainingMandates = remainingMandates + 1;
                var republicVotingNumber = (totalRemainders + fixedRemainingMandates / 2) / fixedRemainingMandates;

                // Na tomto základě se přikáže každé politické straně, politickému hnutí a koalici tolik mandátů, kolikrát je republikové volební číslo obsaženo v součtu zbytků hlasů odevzdaných pro
                // jednotlivou politickou stranu, politické hnutí nebo koalici.
                var partyDivision = partyRemainders.ToDictionary(pkv => pkv.Key, pkv => (quot: Math.DivRem(pkv.Value, republicVotingNumber, out var rem), rem));
                foreach (var pkv in partyDivision)
                {
                    partyMandates.TryGetValue(pkv.Key, out var currMandates);
                    partyMandates[pkv.Key] = currMandates + pkv.Value.quot;
                    remainingMandates -= pkv.Value.quot;
                }

                if (remainingMandates > 0)
                {
                    // (3) Nebyly-li tímto způsobem přikázány všechny mandáty, přikáží se zbývající mandáty postupně těm politickým stranám, politickým hnutím a koalicím, které vykazují největší zbytek dělení
                    // podle odstavce 2;
                    // TODO: Stejně se postupuje, jestliže politická strana, politické hnutí nebo koalice má pro druhé skrutinium méně kandidátů, než kolik mandátů na ni připadá.

                    foreach (var luckyPartyRemGroup in partyDivision.GroupBy(pd => pd.Value.rem).OrderByDescending(g => g.Key))
                    {
                        // při rovnosti zbytků se přikáže mandát té politické straně, politickému hnutí nebo koalici, která má větší součet zbytků hlasů převáděných do druhého skrutinia.
                        foreach (var luckyPartyTotalRemGroup in luckyPartyRemGroup.GroupBy(lprg => partyRemainders[lprg.Key]).OrderByDescending(g => g.Key))
                        {
                            // Jsou-li tyto součty zbytků hlasů stejné, přikáže se mandát té politické straně, politickému hnutí nebo koalici, která obdržela větší počet hlasů;
                            foreach (var luckyPartyTotalVotesGroup in luckyPartyTotalRemGroup.GroupBy(lptvg => partyTotalVotes![lptvg.Key]).OrderByDescending(g => g.Key))
                            {
                                var group = luckyPartyTotalVotesGroup.ToList();
                                // je-li tento počet hlasů stejný, rozhodne los.
                                var chosen = rng.Next(group.Count);
                                var luckyParty = group[chosen].Key;
                                ++partyMandates[luckyParty];
                                --remainingMandates;

                                if (remainingMandates == 0) break;
                            }
                            if (remainingMandates == 0) break;
                        }
                        if (remainingMandates == 0) break;
                    }
                    Debug.Assert(remainingMandates == 0);
                }

                // (4) Bylo-li takto přikázáno o jeden mandát více, než se mělo podle odstavce 2 přikázat,
                if (remainingMandates < 0)
                {
                    Debug.Assert(remainingMandates == -1);

                    // odečte se přebývající mandát té politické straně, politickému hnutí nebo koalici, která vykázala nejmenší zbytek dělení ve druhém skrutiniu.
                    var unluckyPartyRemGroup = partyDivision.GroupBy(pd => pd.Value.rem).OrderBy(g => g.Key).First();
                    // Při stejném zbytku se odečte přebývající mandát té politické straně, politickému hnutí nebo koalici, která získala menší počet hlasů;
                    var unluckyPartyTotalVotesGroup = unluckyPartyRemGroup.GroupBy(ulptvg => partyTotalVotes![ulptvg.Key]).OrderByDescending(g => g.Key).First();
                    var group = unluckyPartyTotalVotesGroup.ToList();
                    // je-li tento počet hlasů stejný, rozhodne los.                
                    var chosen = rng.Next(group.Count);
                    var unluckyParty = group[chosen].Key;
                    --partyMandates[unluckyParty];
                    ++remainingMandates;
                    Debug.Assert(partyMandates[unluckyParty] >= 0);
                }
            }

            Debug.Assert(remainingMandates == 0);
            return partyMandates;
        }

        private static Dictionary<Region, Dictionary<Party, int>> GenerateRandomResults(Dictionary<Region, int> citizens, Dictionary<Region, RandomVariableCharacteristics> expectedTurnout, Dictionary<Region, Dictionary<Party, RandomVariableCharacteristics>> expectedResults) =>
            expectedResults
                .ToDictionary(rkv => rkv.Key,
                    rkv =>
                    {
                        var region = rkv.Key;
                        var turnout = GenerateBinomialRandomNumber(citizens[region], expectedTurnout[region].Mean, expectedTurnout[region].StdDev * 10);
                        Debug.Assert(turnout > 0);
                        // var generatedResults = expectedResults[region].ToDictionary(pkv => pkv.Key, pkv => turnout * GenerateLogNormalRandomNumber(pkv.Value));
                        var generatedResults = expectedResults[region].ToDictionary(pkv => pkv.Key, pkv => GenerateBinomialRandomNumber(turnout, pkv.Value.Mean, pkv.Value.StdDev * 10));
                        Debug.Assert(generatedResults.Values.All(r => r >= 0));
                        var generatedTotal = generatedResults.Values.Sum();
                        return generatedResults.ToDictionary(pkv => pkv.Key, pkv => (int) Math.Round((double) pkv.Value * turnout / generatedTotal));

                        // var turnout = MathNet.Numerics.Distributions.Binomial.Sample(rng, expectedTurnout[region].Mean, citizens[region]);
                        // var generatedResults = MathNet.Numerics.Distributions.Multinomial.Sample(rng, expectedResults[region].Values.Select(v => v.Mean).ToArray(), turnout);
                        // return expectedResults[region].Select((pkv, i) => (party: pkv.Key, result: generatedResults[i])).ToDictionary(e => e.party, e => e.result);
                    });

        private static int GenerateBinomialRandomNumber(int n, double p, double varRate)
        {
            var np = n * p;
            var gaussian = GenerateGaussianRandomNumber(new(np, varRate * Math.Sqrt(np * (1 - p))));
            return Math.Clamp((int) Math.Round(gaussian), 0, n);
        }

        private static double GenerateBinomialCdf(int n, double p, int k) => MathNet.Numerics.Combinatorics.Combinations(n, k) * Math.Pow(p, k) * Math.Pow(1 - p, n - k);

        private static double GenerateLogNormalRandomNumber(RandomVariableCharacteristics chars)
        {
            // var gaussian = GenerateGaussianRandomNumber(new(0, 1));
            // return chars.Mean * Math.Exp(chars.StdDev * gaussian);
            return Math.Exp(GenerateGaussianRandomNumber(new(Math.Log(chars.Mean), chars.StdDev)));
        }

        private static double GenerateGaussianRandomNumber(RandomVariableCharacteristics chars)
        {
            return MathNet.Numerics.Distributions.Normal.Sample(rng, chars.Mean, chars.StdDev);
            // var u1 = 1.0 - rng.NextDouble(); //uniform(0,1] random doubles
            // var u2 = 1.0 - rng.NextDouble();
            // var randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            // return chars.Mean + chars.StdDev * randStdNormal; //random normal(mean,stdDev^2)
        }
    }

    internal class AtomicInteger
    {
        private int value;
        public int Value => value;

        public AtomicInteger(int value)
        {
            this.value = value;
        }

        public void Increment() => Interlocked.Increment(ref value);

        public void Add(int addend) => Interlocked.Add(ref value, addend);

        public void MinWith(int value2)
        {
            while (true)
            {
                var v = value;
                if (value2 >= v) return;

                if (Interlocked.CompareExchange(ref value, value2, v) == v) break;
            }
        }

        public void MaxWith(int value2)
        {
            while (true)
            {
                var v = value;
                if (value2 <= v) return;

                if (Interlocked.CompareExchange(ref value, value2, v) == v) break;
            }
        }

        public override string ToString() => value.ToString();
    }

    public readonly struct RandomVariableCharacteristics
    {
        public readonly double Mean;
        public readonly double StdDev;

        public RandomVariableCharacteristics(double mean, double stdDev)
        {
            Mean = mean;
            StdDev = stdDev;
        }
    }

    public record Region(string Name);

    public record Party(string ShortName);
}