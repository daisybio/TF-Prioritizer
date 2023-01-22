import {NgModule} from '@angular/core';
import {RouterModule, Routes} from '@angular/router';
import {RankingComponent} from "./pages/ranking/ranking.component";
import {ImportantLociComponent} from "./pages/generalAnalysis/important-loci/important-loci.component";
import {TopLog2fcComponent} from "./pages/generalAnalysis/top-log2fc/top-log2fc.component";
import {CoOccurrenceComponent} from "./pages/generalAnalysis/co-occurrence/co-occurrence.component";
import {
  StatisticalEvaluationComponent
} from "./pages/generalAnalysis/statistical-evaluation/statistical-evaluation.component";
import {ParametersComponent} from "./pages/explainatory/parameters/parameters.component";
import {DocumentationComponent} from "./pages/explainatory/documentation/documentation.component";
import {TfDetailsComponent} from "./pages/perTf/tf-details/tf-details.component";
import {
  RegressionCoefficientsComponent
} from "./pages/generalAnalysis/regression-coefficients/regression-coefficients.component";

const routes: Routes = [
  {path: '', component: RankingComponent},
  {path: 'importantLoci', component: ImportantLociComponent},
  {path: 'topLog2fc', component: TopLog2fcComponent},
  {path: 'coOccurrence', component: CoOccurrenceComponent},
  {path: 'statisticalEvaluation', component: StatisticalEvaluationComponent},
  {path: 'importantLoci', component: ImportantLociComponent},
  {path: 'parameters', component: ParametersComponent},
  {path: 'documentation', component: DocumentationComponent},
  {path: 'details', component: TfDetailsComponent},
  {path: 'regressionCoefficients', component: RegressionCoefficientsComponent},

  {path: '**', redirectTo: ''}
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule {
}
