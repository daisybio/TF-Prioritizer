import {Routes, RouterModule} from '@angular/router';

import {TfRankingComponent} from "./tf-ranking/tf-ranking.component";
import {ValidationComponent} from "./validation/validation.component";
import {NgModule} from "@angular/core";

const routes: Routes = [
  {path: '', component: TfRankingComponent},
  {path: 'validation', component: ValidationComponent},

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
