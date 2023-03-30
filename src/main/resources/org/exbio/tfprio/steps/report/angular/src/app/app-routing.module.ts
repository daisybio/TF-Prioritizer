import {NgModule} from '@angular/core';
import {RouterModule, Routes} from '@angular/router';
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
import {HomeComponent} from "./pages/home/home.component";

const routes: Routes = [
  {path: '', component: HomeComponent},
  {path: 'coOccurrence', component: CoOccurrenceComponent},
  {path: 'statisticalEvaluation', component: StatisticalEvaluationComponent},
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
