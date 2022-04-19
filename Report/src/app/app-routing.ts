import {Routes, RouterModule} from '@angular/router';

import {TfRankingComponent} from "./tf-ranking/tf-ranking.component";
import {ValidationComponent} from "./validation/validation.component";
import {NgModule} from "@angular/core";
import {RegressionComponent} from "./regression/regression.component";
import {DistributionComponent} from "./distribution/distribution.component";
import {DocumentationComponent} from "./documentation/documentation.component";
import {ParametersComponent} from "./parameters/parameters.component";

const routes: Routes = [
  {path: '', component: TfRankingComponent},
  {path: 'validation', component: ValidationComponent},
  {path: 'regression', component: RegressionComponent},
  {path: 'distribution', component: DistributionComponent},
  {path: 'documentation', component: DocumentationComponent},
  {path: 'parameters', component: ParametersComponent},

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
