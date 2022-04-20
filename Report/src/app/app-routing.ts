import {Routes, RouterModule} from '@angular/router';

import {TfRankingComponent} from "./tf-ranking/tf-ranking.component";
import {ValidationComponent} from "./validation/validation.component";
import {NgModule} from "@angular/core";
import {RegressionComponent} from "./regression/regression.component";
import {DistributionComponent} from "./distribution/distribution.component";
import {DocumentationComponent} from "./documentation/documentation.component";
import {ParametersComponent} from "./parameters/parameters.component";
import {ImportantLociComponent} from "./important-loci/important-loci.component";
import {TopLog2fcComponent} from "./top-log2fc/top-log2fc.component";
import {CoOccurrenceAnalysisComponent} from "./co-occurrence-analysis/co-occurrence-analysis.component";

const routes: Routes = [
  {path: '', component: TfRankingComponent},
  {path: 'validation', component: ValidationComponent},
  {path: 'regression', component: RegressionComponent},
  {path: 'distribution', component: DistributionComponent},
  {path: 'documentation', component: DocumentationComponent},
  {path: 'parameters', component: ParametersComponent},
  {path: 'importantLoci', component: ImportantLociComponent},
  {path: 'topLog2fc', component: TopLog2fcComponent},
  {path: 'coOccurrenceAnalysis', component: CoOccurrenceAnalysisComponent},

  // otherwise redirect to home
  {path: '**', redirectTo: ''}
];

// configures NgModule imports and exports
@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule {
}
